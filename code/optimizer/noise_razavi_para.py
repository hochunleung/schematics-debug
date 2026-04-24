import numpy as np
import pandas as pd
import re
from lut_engine import MosData
from scipy.interpolate import RegularGridInterpolator

# ================= 核心求解流程 =================
import numpy as np

def evaluate_differential_noise(op_config, freqs, T=300):
    """
    全差分運放 (Buffer 配置) 雜訊評估函數 - 子矩陣精確求解版
    :param op_config: 操作點字典，包含 'model' (如 nch, pch 等 LUT 實例)
    :param freqs: 頻率陣列或列表
    :return: 包含頻率、增益、總輸出雜訊、等效輸入雜訊及各元件貢獻的字典
    """
    # ==========================================
    # 1. 定義節點 (擴充外部理想輸入節點)
    # ==========================================
    N_P1, N_N1, N_S = 0, 1, 2
    N_1, N_2, N_3 = 3, 4, 5
    N_P2, N_N2, N_4 = 6, 7, 8
    # 新增：理想電壓源節點 (Inputs)
    N_IN_P, N_IN_N = 9, 10 
    
    NUM_NODES = 11
    free_nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8]
    inputs = [9, 10]

    # ==========================================
    # 2. 提取參數
    # ==========================================
    K_dict, C_dict = {}, {}
    for inst, cfg in op_config['mosfets'].items():
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
    kT4 = 4 * 1.380649e-23 * T

    # 準備儲存結果的容器
    results = {
        'freqs': freqs, 'gain': [], 'out_noise_psd': [], 'irn_psd': [],
        'contributions': {comp: [] for comp in list(op_config['mosfets'].keys()) + ['R_in', 'R_fb', 'R_cmfb', 'R_miller']}
    }
    results['contributions'].pop('passives', None)

    # ==========================================
    # 3. 頻率迭代與矩陣求解
    # ==========================================
    for f in freqs:
        w = 2 * np.pi * f
        s = 1j * w if w != 0 else 1e-15j
        
        y_miller = 1.0 / (passives['Rz'] + 1.0/(s * passives['Cc']))
        
        # --- A. 建立 Y 矩陣 ---
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

        # 被動網路 (精確分離 Input 與 Feedback)
        R_in, R_cmfb = passives['R_in'], passives['R_cmfb']
        # 1. 輸入電阻 (連接理想源與 Gate)
        add_y(N_P1, N_IN_P, 1/R_in); add_y(N_N1, N_IN_N, 1/R_in)
        # 2. 回饋電阻 (連接 Gate 與 輸出)
        add_y(N_P1, N_N2, 1/R_in); add_y(N_N1, N_P2, 1/R_in)
        
        add_y(N_P2, N_4, 1/R_cmfb); add_y(N_N2, N_4, 1/R_cmfb)
        add_y(N_1, N_3, 1/R_cmfb); add_y(N_2, N_3, 1/R_cmfb)
        add_y(N_P2, N_1, y_miller); add_y(N_N2, N_2, y_miller)

        # 寄生電容與主動元件 (與原本相同)
        add_y(N_P1, N_S, s*C_dict['M1']['Cgs']); add_y(N_P1, N_1, s*C_dict['M1']['Cgd'])
        add_y(N_N1, N_S, s*C_dict['M2']['Cgs']); add_y(N_N1, N_2, s*C_dict['M2']['Cgd'])
        add_y(N_3, -1, s*C_dict['M3']['Cgs']); add_y(N_3, N_1, s*C_dict['M3']['Cgd'])
        add_y(N_3, -1, s*C_dict['M4']['Cgs']); add_y(N_3, N_2, s*C_dict['M4']['Cgd'])
        add_y(N_1, -1, s*C_dict['M5']['Cgs']); add_y(N_1, N_P2, s*C_dict['M5']['Cgd']); add_y(N_P2, -1, s*C_dict['M5']['Cdb'])
        add_y(N_2, -1, s*C_dict['M6']['Cgs']); add_y(N_2, N_N2, s*C_dict['M6']['Cgd']); add_y(N_N2, -1, s*C_dict['M6']['Cdb'])
        add_y(N_4, -1, s*C_dict['M7']['Cgs']); add_y(N_4, N_P2, s*C_dict['M7']['Cgd']); add_y(N_P2, -1, s*C_dict['M7']['Cdb'])
        add_y(N_4, -1, s*C_dict['M8']['Cgs']); add_y(N_4, N_N2, s*C_dict['M8']['Cgd']); add_y(N_N2, -1, s*C_dict['M8']['Cdb'])
        add_y(N_S, -1, s*(C_dict['M10']['Cdb'] + C_WELL_S))

        add_gm(N_1, N_P1, N_S, K_dict['M1']['K10']); add_y(N_1, N_S, K_dict['M1']['K01'])
        add_gm(N_2, N_N1, N_S, K_dict['M2']['K10']); add_y(N_2, N_S, K_dict['M2']['K01'])
        add_gm(N_1, N_3, -1, K_dict['M3']['K10']);   add_y(N_1, -1, K_dict['M3']['K01'])
        add_gm(N_2, N_3, -1, K_dict['M4']['K10']);   add_y(N_2, -1, K_dict['M4']['K01'])
        add_gm(N_P2, N_1, -1, K_dict['M5']['K10']);  add_y(N_P2, -1, K_dict['M5']['K01'])
        add_gm(N_N2, N_2, -1, K_dict['M6']['K10']);  add_y(N_N2, -1, K_dict['M6']['K01'])
        add_gm(N_P2, N_4, -1, K_dict['M7']['K10']);  add_y(N_P2, -1, K_dict['M7']['K01'])
        add_gm(N_N2, N_4, -1, K_dict['M8']['K10']);  add_y(N_N2, -1, K_dict['M8']['K01'])
        add_y(N_S, -1, K_dict['M10']['K01'])

        # --- B. 建立 SI 獨立貢獻字典 ---
        SI_dict = {comp: np.zeros((NUM_NODES, NUM_NODES)) for comp in results['contributions'].keys()}
        
        def add_i2(mat, n1, n2, i2):
            if n1 >= 0: mat[n1, n1] += i2
            if n2 >= 0: mat[n2, n2] += i2
            if n1 >= 0 and n2 >= 0:
                mat[n1, n2] -= i2; mat[n2, n1] -= i2

        # 1. 填寫 MOSFET 雜訊
        for inst, cfg in op_config['mosfets'].items():
            if inst == 'passives': continue
            m, L, W, VGS, VDS = cfg['model'], cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
            i2_n = m.get_noise_psd(L, VGS, VDS, W, f)
            
            mat = SI_dict[inst]
            if inst == 'M1': add_i2(mat, N_1, N_S, i2_n)
            elif inst == 'M2': add_i2(mat, N_2, N_S, i2_n)
            elif inst == 'M3': add_i2(mat, N_1, -1, i2_n)
            elif inst == 'M4': add_i2(mat, N_2, -1, i2_n)
            elif inst == 'M5': add_i2(mat, N_P2, -1, i2_n)
            elif inst == 'M6': add_i2(mat, N_N2, -1, i2_n)
            elif inst == 'M7': add_i2(mat, N_P2, -1, i2_n)
            elif inst == 'M8': add_i2(mat, N_N2, -1, i2_n)
            elif inst == 'M10': add_i2(mat, N_S, -1, i2_n)

        # 2. 填寫被動元件雜訊 (分離 R_in 與 R_fb)
        add_i2(SI_dict['R_in'], N_P1, N_IN_P, kT4 / R_in)
        add_i2(SI_dict['R_in'], N_N1, N_IN_N, kT4 / R_in)
        
        add_i2(SI_dict['R_fb'], N_P1, N_N2, kT4 / R_in)
        add_i2(SI_dict['R_fb'], N_N1, N_P2, kT4 / R_in)
        
        add_i2(SI_dict['R_cmfb'], N_P2, N_4, kT4 / R_cmfb)
        add_i2(SI_dict['R_cmfb'], N_N2, N_4, kT4 / R_cmfb)
        add_i2(SI_dict['R_cmfb'], N_1, N_3, kT4 / R_cmfb)
        add_i2(SI_dict['R_cmfb'], N_2, N_3, kT4 / R_cmfb)

        i2_miller = kT4 * passives['Rz'] * (np.abs(y_miller)**2)
        add_i2(SI_dict['R_miller'], N_P2, N_1, i2_miller)
        add_i2(SI_dict['R_miller'], N_N2, N_2, i2_miller)

        # --- C. 提取子矩陣與求解 ---
        Y_ff = Y[np.ix_(free_nodes, free_nodes)]
        Y_fi = Y[np.ix_(free_nodes, inputs)]
        
        try:
            Z_ff = np.linalg.inv(Y_ff)
        except np.linalg.LinAlgError:
            continue

        # 1. 求解閉環增益 (設定差分輸入 V_IN_P = 0.5V, V_IN_N = -0.5V)
        # V_in = np.array([0.5, -0.5]) 
        # V_free = -Z_ff @ (Y_fi @ V_in)
        # 1. 求解閉環增益 (Robust 版本)
        # 建立「全節點規模」的輸入電壓向量
        V_in_full = np.zeros(NUM_NODES)
        V_in_full[N_IN_P] = 0.5
        V_in_full[N_IN_N] = -0.5
        
        # 自動提取屬於 inputs 節點的激勵電壓
        V_in = V_in_full[inputs] 
        
        V_free = -Z_ff @ (Y_fi @ V_in)
        
        # 利用同樣的切片邏輯計算差分輸出 (避免 V_free 索引錯位)
        # 建立全節點的電壓分佈 (已知 inputs 節點電壓，結合求出的 free 節點電壓)
        V_nodes_full = np.zeros(NUM_NODES, dtype=complex)
        V_nodes_full[free_nodes] = V_free
        V_nodes_full[inputs] = V_in
        # 差分輸出電壓 = V(P2) - V(N2)
        diff_gain = np.abs(V_free[N_P2] - V_free[N_N2])
        results['gain'].append(diff_gain)

        # 2. 求解各元件的雜訊貢獻
        # 建立差分輸出提取向量 [0,0,..., 1(P2), -1(N2), 0]
        # c_out = np.zeros(len(free_nodes))
        # c_out[N_P2] = 1; c_out[N_N2] = -1
        # ==========================================
        # 2. 求解各元件的雜訊貢獻 (Robust 版本)
        # ==========================================
        # 建立「全節點規模」的差分輸出提取向量
        c_out_full = np.zeros(NUM_NODES)
        c_out_full[N_P2] = 1
        c_out_full[N_N2] = -1
        
        # 透過 Numpy 陣列切片，自動精準映射到 free_nodes 的相對維度！
        c_out = c_out_full[free_nodes]
        # H_out 是將所有 free_nodes 電流轉換到差分輸出電壓的轉移向量
        H_out = c_out @ Z_ff 

        total_out_psd = 0
        for comp_name, SI_mat in SI_dict.items():
            SI_ff = SI_mat[np.ix_(free_nodes, free_nodes)]
            # 向量-矩陣-向量 乘法得到純量 PSD
            comp_out_psd = np.real(H_out @ SI_ff @ H_out.conj().T)
            results['contributions'][comp_name].append(comp_out_psd)
            total_out_psd += comp_out_psd
            
        results['out_noise_psd'].append(total_out_psd)
        
        # 3. 計算等效輸入雜訊 (IRN)
        results['irn_psd'].append(total_out_psd / (diff_gain**2))

    return results

if __name__ == "__main__":
    # 載入 LUT 資料庫 (注意 P 管開啟 is_pmos=True)
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
    op_config = {
        'mosfets': {
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
        },
        'passives': {'R_f': 4000, 'R_in': 4000, 'Rz': 204.6, 'Cc': 402e-15, 'R_cmfb': 40e3}
    }
    noise = evaluate_differential_noise(op_config, [10e6])    
    print(f"Gain      : {noise['gain'][0]} ")
    print(f"out noise : {noise['out_noise_psd'][0]} ")
    print(f"in  noise : {noise['irn_psd'][0]} ")
    