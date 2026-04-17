import numpy as np
import pandas as pd
import re
from lut_engine import MosData
from scipy.interpolate import RegularGridInterpolator

# ================= 核心求解流程 =================
import numpy as np

def evaluate_differential_noise(op_config, freqs, T=300):
    """
    全差分運放 Noise 評估核心函數
    :param op_config: 來自 map_vars_to_params 的操作點字典，包含 'model' 實例
    :param freqs: 頻率陣列 (Numpy array)
    :param T: 絕對溫度 (預設 300K)
    :return: dict 包含頻率、差分輸出雜訊 PSD、差分輸入等效雜訊 PSD 與閉迴路增益
    """
    # 節點索引定義
    N_P1, N_N1, N_S = 0, 1, 2
    N_1, N_2, N_3 = 3, 4, 5
    N_P2, N_N2, N_4 = 6, 7, 8
    NUM_NODES = 9

    # 1. 動態提取 K 字典和寄生電容字典 (與原邏輯相同)
    K_dict = {}
    C_dict = {}
    for inst, cfg in op_config.items():
        if inst == 'passives': continue
        m = cfg['model']
        L, W, VGS, VDS = cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
        
        K_dict[inst] = {k: m.lookup_by_vgs(k+'_W', L, VGS, VDS) * W for k in ['K10', 'K01']}
        C_dict[inst] = {
            'Cgs': m.lookup_by_vgs('CGS_W', L, VGS, VDS) * W,
            'Cgd': m.lookup_by_vgs('CGD_W', L, VGS, VDS) * W,
            'Cdb': np.abs(m.lookup_by_vgs('CDB_W', L, VGS, VDS) * W) 
        }

    C_WELL_S = 3 * C_dict['M10']['Cdb'] 
    passives = op_config['passives']
    kT4 = 4 * 1.380649e-23 * T  # 熱雜訊常數 4kT

    # 2. 構建動態 Y 矩陣 (保留你的原始結構，僅由外部傳入 s 與 y_miller)
    def build_Y_matrix(s, y_miller):
        Y = np.zeros((NUM_NODES, NUM_NODES), dtype=complex)
        
        def add_y(n1, n2, y_val):
            if n1 >= 0: Y[n1, n1] += y_val
            if n2 >= 0: Y[n2, n2] += y_val
            if n1 >= 0 and n2 >= 0:
                Y[n1, n2] -= y_val; Y[n2, n1] -= y_val

        def add_gm(d, g, src, gm):
            if d >= 0:
                if g >= 0: Y[d, g] += gm
                if src >= 0: Y[d, src] -= gm
            if src >= 0:
                if g >= 0: Y[src, g] -= gm
                if src >= 0: Y[src, src] += gm

        # --- 被動元件 ---
        R_in, R_cmfb = passives['R_in'], passives['R_cmfb']
        add_y(N_P1, N_N2, 1/R_in); add_y(N_N1, N_P2, 1/R_in)
        add_y(N_P1, -1, 1/R_in); add_y(N_N1, -1, 1/R_in)
        add_y(N_P2, N_4, 1/R_cmfb); add_y(N_N2, N_4, 1/R_cmfb)
        add_y(N_1, N_3, 1/R_cmfb); add_y(N_2, N_3, 1/R_cmfb)
        add_y(N_P2, N_1, y_miller); add_y(N_N2, N_2, y_miller)

        # --- 電容注入 ---
        add_y(N_P1, N_S, s*C_dict['M1']['Cgs']); add_y(N_P1, N_1, s*C_dict['M1']['Cgd'])
        add_y(N_N1, N_S, s*C_dict['M2']['Cgs']); add_y(N_N1, N_2, s*C_dict['M2']['Cgd'])
        add_y(N_3, -1, s*C_dict['M3']['Cgs']); add_y(N_3, N_1, s*C_dict['M3']['Cgd'])
        add_y(N_3, -1, s*C_dict['M4']['Cgs']); add_y(N_3, N_2, s*C_dict['M4']['Cgd'])
        add_y(N_1, -1, s*C_dict['M5']['Cgs']); add_y(N_1, N_P2, s*C_dict['M5']['Cgd']); add_y(N_P2, -1, s*C_dict['M5']['Cdb'])
        add_y(N_2, -1, s*C_dict['M6']['Cgs']); add_y(N_2, N_N2, s*C_dict['M6']['Cgd']); add_y(N_N2, -1, s*C_dict['M6']['Cdb'])
        add_y(N_4, -1, s*C_dict['M7']['Cgs']); add_y(N_4, N_P2, s*C_dict['M7']['Cgd']); add_y(N_P2, -1, s*C_dict['M7']['Cdb'])
        add_y(N_4, -1, s*C_dict['M8']['Cgs']); add_y(N_4, N_N2, s*C_dict['M8']['Cgd']); add_y(N_N2, -1, s*C_dict['M8']['Cdb'])
        add_y(N_S, -1, s*(C_dict['M10']['Cdb'])) # + C_WELL_S

        # --- 跨導與輸出電導 ---
        add_gm(N_1, N_P1, N_S, K_dict['M1']['K10']); add_y(N_1, N_S, K_dict['M1']['K01'])
        add_gm(N_2, N_N1, N_S, K_dict['M2']['K10']); add_y(N_2, N_S, K_dict['M2']['K01'])
        add_gm(N_1, N_3, -1, K_dict['M3']['K10']);   add_y(N_1, -1, K_dict['M3']['K01'])
        add_gm(N_2, N_3, -1, K_dict['M4']['K10']);   add_y(N_2, -1, K_dict['M4']['K01'])
        add_gm(N_P2, N_1, -1, K_dict['M5']['K10']);  add_y(N_P2, -1, K_dict['M5']['K01'])
        add_gm(N_N2, N_2, -1, K_dict['M6']['K10']);  add_y(N_N2, -1, K_dict['M6']['K01'])
        add_gm(N_P2, N_4, -1, K_dict['M7']['K10']);  add_y(N_P2, -1, K_dict['M7']['K01'])
        add_gm(N_N2, N_4, -1, K_dict['M8']['K10']);  add_y(N_N2, -1, K_dict['M8']['K01'])
        add_y(N_S, -1, K_dict['M10']['K01'])

        return Y

    # 3. 構建動態雜訊電流矩陣 SI
    def build_SI_matrix(f, y_miller):
        SI = np.zeros((NUM_NODES, NUM_NODES))
        
        def add_i2(n1, n2, i2):
            if n1 >= 0: SI[n1, n1] += i2
            if n2 >= 0: SI[n2, n2] += i2
            if n1 >= 0 and n2 >= 0:
                SI[n1, n2] -= i2; SI[n2, n1] -= i2

        # --- A. MOSFET 雜訊注入 (呼叫 LUT，並對應你的 Drain / Source 節點) ---
        for inst, cfg in op_config.items():
            if inst == 'passives': continue
            m = cfg['model']
            L, W, VGS, VDS = cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
            
            # 從 LUT 取得該頻率的電流 PSD
            i2_n = m.get_noise_psd(L, VGS, VDS, W, f)
            
            if inst == 'M1': add_i2(N_1, N_S, i2_n)
            elif inst == 'M2': add_i2(N_2, N_S, i2_n)
            elif inst == 'M3': add_i2(N_1, -1, i2_n)
            elif inst == 'M4': add_i2(N_2, -1, i2_n)
            elif inst == 'M5': add_i2(N_P2, -1, i2_n)
            elif inst == 'M6': add_i2(N_N2, -1, i2_n)
            elif inst == 'M7': add_i2(N_P2, -1, i2_n)
            elif inst == 'M8': add_i2(N_N2, -1, i2_n)
            elif inst == 'M10': add_i2(N_S, -1, i2_n)

        # --- B. 被動元件熱雜訊注入 ---
        R_in, R_cmfb = passives['R_in'], passives['R_cmfb']
        add_i2(N_P1, N_N2, kT4 / R_in); add_i2(N_N1, N_P2, kT4 / R_in)
        add_i2(N_P1, -1, kT4 / R_in); add_i2(N_N1, -1, kT4 / R_in)
        
        add_i2(N_P2, N_4, kT4 / R_cmfb); add_i2(N_N2, N_4, kT4 / R_cmfb)
        add_i2(N_1, N_3, kT4 / R_cmfb); add_i2(N_2, N_3, kT4 / R_cmfb)

        # --- C. Miller 電阻的等效 Norton 雜訊電流 ---
        # 串聯 R-C 分支的等效電流雜訊為 I_n^2 = 4kT * Rz * |Y_miller|^2
        i2_miller = kT4 * passives['Rz'] * (np.abs(y_miller)**2)
        add_i2(N_P2, N_1, i2_miller)
        add_i2(N_N2, N_2, i2_miller)

        return SI

    # 4. 準備求解矩陣與儲存結果
    S_out_diff = np.zeros(len(freqs))
    S_in_diff = np.zeros(len(freqs))
    Av_mag = np.zeros(len(freqs))

    # 定義差分提取向量
    c_out = np.zeros(NUM_NODES); c_out[N_P2] = 1; c_out[N_N2] = -1
    c_in = np.zeros(NUM_NODES); c_in[N_P1] = 1; c_in[N_N1] = -1

    for i, f in enumerate(freqs):
        w = 2 * np.pi * f
        s = 1j * w if w != 0 else 1e-15j # 防止 DC 除以 0
        
        y_miller = 1.0 / (passives['Rz'] + 1.0/(s * passives['Cc']))
        
        # 建立導納與雜訊矩陣
        Y = build_Y_matrix(s, y_miller)
        SI = build_SI_matrix(f, y_miller)
        
        # 求解阻抗矩陣 Z = Y^-1
        Z = np.linalg.inv(Y)
        
        # --- 步驟 A：計算差分輸出雜訊 PSD ---
        Z_out_diff = c_out @ Z  # 提取輸出差分行的阻抗向量
        s_out = (Z_out_diff @ SI @ Z_out_diff.conj().T).real
        S_out_diff[i] = s_out
        
        # --- 步驟 B：利用 MNA 推導精確的閉迴路電壓增益 ---
        # 藉由在差分輸入端注入 1A 測試電流，反推各節點電壓
        V_nodes_test = Z @ c_in.T  
        V_in_diff_test = c_in @ V_nodes_test   # 考慮了 Rin 與 Cgs 後的實際輸入跨壓
        V_out_diff_test = c_out @ V_nodes_test # 實際產生的輸出電壓
        
        Av = V_out_diff_test / V_in_diff_test
        Av_mag[i] = np.abs(Av)
        
        # --- 步驟 C：計算輸入等效雜訊 PSD ---
        S_in_diff[i] = s_out / (np.abs(Av)**2)

    return {
        'freqs': freqs,
        'S_out_diff': S_out_diff,
        'S_in_diff': S_in_diff,
        'Gain_V_V': Av_mag
    }


if __name__ == "__main__":
    # ⚠️ 請填入該全差分 Op-amp 在 Cadence dcOp 中的實際電壓 (PMOS 請填絕對值)
    # OP_CONFIG = {
    #     'M1':  {'vgs': 0.2565, 'vds': 0.4191, 'model': 'nch_lvt', 'type_inl': 'nch'}, # 第一級輸入對管 (左)
    #     'M2':  {'vgs': 0.2565, 'vds': 0.4191, 'model': 'nch_lvt', 'type_inl': 'nch'}, # 第一級輸入對管 (右)
    #     'M3':  {'vgs': 0.2873, 'vds': 0.2873, 'model': 'pch', 'type_inl': 'pch'},     # 第一級負載 (左)
    #     'M4':  {'vgs': 0.2873, 'vds': 0.2873, 'model': 'pch', 'type_inl': 'pch'},     # 第一級負載 (右)
    #     'M5':  {'vgs': 0.2873, 'vds': 0.4498, 'model': 'pch', 'type_inl': 'pch'},    # 第二級驅動 (左)
    #     'M6':  {'vgs': 0.2873, 'vds': 0.4498, 'model': 'pch', 'type_inl': 'pch'},    # 第二級驅動 (右)
    #     'M7':  {'vgs': 0.3301, 'vds': 0.4502, 'model': 'nch', 'type_inl': 'nch'},     # 第二級負載 (左)
    #     'M8':  {'vgs': 0.3301, 'vds': 0.4502, 'model': 'nch', 'type_inl': 'nch'},     # 第二級負載 (右)
    #     'M10': {'vgs': 0.4585, 'vds': 0.1936, 'model': 'nch', 'type_inl': 'nch'}     # 尾電流源
    # }
    
    # 載入 LUT 資料庫 (注意 P 管開啟 is_pmos=True)
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
    op_config = {
        'M1':  {'VGS': 0.2565, 'VDS': 0.4191, 'W': 96e-6, 'L': 100e-9, 'model': nch_lvt, 'type_inl': 'nch'},
        'M2':  {'VGS': 0.2565, 'VDS': 0.4191, 'W': 96e-6, 'L': 100e-9, 'model': nch_lvt, 'type_inl': 'nch'},
        'M3':  {'VGS': 0.2873, 'VDS': 0.2873, 'W': 96e-6, 'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
        'M4':  {'VGS': 0.2873, 'VDS': 0.2873, 'W': 96e-6, 'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
        'M5':  {'VGS': 0.2873, 'VDS': 0.4498, 'W': 192e-6,'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
        'M6':  {'VGS': 0.2873, 'VDS': 0.4498, 'W': 192e-6,'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
        'M7':  {'VGS': 0.3301, 'VDS': 0.4502, 'W': 48e-6, 'L': 100e-9, 'model': nch, 'type_inl': 'nch'},
        'M8':  {'VGS': 0.3301, 'VDS': 0.4502, 'W': 48e-6, 'L': 100e-9, 'model': nch, 'type_inl': 'nch'},
        'M10': {'VGS': 0.4585, 'VDS': 0.1936, 'W': 32e-6, 'L': 500e-9, 'model': nch, 'type_inl': 'nch'},
        # 被動元件一併傳遞
        'passives': {'R_f': 4000, 'R_in': 4000, 'Rz': 204.6, 'Cc': 402e-15, 'R_cmfb': 40e3}
    }
    noise = evaluate_differential_noise(op_config, [10e3])
    # metrics = evaluate_volterra_sfdr(op_config, V_AMP_IN=0.316, FIN=53/512*100e6)
    
    print(f"out noise : {noise['S_out_diff'][0]} ")