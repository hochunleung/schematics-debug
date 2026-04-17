import numpy as np
import pandas as pd
import re
from lut_engine import MosData
from scipy.interpolate import RegularGridInterpolator

# ================= 核心求解流程 =================
def evaluate_volterra_sfdr(op_config, V_AMP_IN=0.632, FIN=10.35e6):
    """
    Volterra 優化器核心函數
    :param op_config: 來自 map_vars_to_params 的操作點字典
    :return: dict 包含 A1, HD2_dB, HD3_dB, SFDR
    """
    # 節點索引定義
    N_P1, N_N1, N_S = 0, 1, 2
    N_1, N_2, N_3 = 3, 4, 5
    N_P2, N_N2, N_4 = 6, 7, 8

    # 1. 動態提取 K 字典和寄生電容字典
    K_dict = {}
    C_dict = {}
    for inst, cfg in op_config.items():
        if inst == 'passives': continue
        m = cfg['model']
        L, W, VGS, VDS = cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
        
        # 提取 K 參數 (密度 * 實際寬度)
        K_dict[inst] = {k: m.lookup_by_vgs(k+'_W', L, VGS, VDS) * W for k in ['K10', 'K01', 'K20', 'K02', 'K11', 'K30', 'K03', 'K21', 'K12']}
        print(f"[{inst}] gm: {K_dict[inst]['K10']*1000:.3f} mS | gds: {K_dict[inst]['K01']*1e6:.2f} uS")
        
        # 提取寄生電容 (密度 * 實際寬度，並確保 Cdb 為正絕對值)
        C_dict[inst] = {
            'Cgs': m.lookup_by_vgs('CGS_W', L, VGS, VDS) * W,
            'Cgd': m.lookup_by_vgs('CGD_W', L, VGS, VDS) * W,
            'Cdb': np.abs(m.lookup_by_vgs('CDB_W', L, VGS, VDS) * W) 
        }

    # 預估井電容 (例如以 M10 的 Cdb 按比例估算，可視需要調整)
    C_WELL_S = 3 * C_dict['M10']['Cdb'] 
    
    passives = op_config['passives']

    # 2. 構建動態 Y 矩陣
    def build_Y_matrix(w):
        Y = np.zeros((9, 9), dtype=complex)
        s = 1j * w
        
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

        # 被動元件
        R_in, R_cmfb = passives['R_in'], passives['R_cmfb']
        add_y(N_P1, N_N2, 1/R_in); add_y(N_N1, N_P2, 1/R_in)
        add_y(N_P1, -1, 1/R_in); add_y(N_N1, -1, 1/R_in)
        add_y(N_P2, N_4, 1/R_cmfb); add_y(N_N2, N_4, 1/R_cmfb)
        add_y(N_1, N_3, 1/R_cmfb); add_y(N_2, N_3, 1/R_cmfb)
        
        y_miller = 1.0 / (passives['Rz'] + 1.0/(s * passives['Cc']))
        add_y(N_P2, N_1, y_miller); add_y(N_N2, N_2, y_miller)

        # 動態電容注入 (包含 S-B 短接拓撲)
        # M1/M2 (Bulk=Source -> Cdb 變成 Cds 跨接)
        add_y(N_P1, N_S, s*C_dict['M1']['Cgs']); add_y(N_P1, N_1, s*C_dict['M1']['Cgd'])#; add_y(N_1, N_S, s*C_dict['M1']['Cdb'])
        add_y(N_N1, N_S, s*C_dict['M2']['Cgs']); add_y(N_N1, N_2, s*C_dict['M2']['Cgd'])#; add_y(N_2, N_S, s*C_dict['M2']['Cdb'])
        
        # M3/M4
        add_y(N_3, -1, s*C_dict['M3']['Cgs']); add_y(N_3, N_1, s*C_dict['M3']['Cgd'])#; add_y(N_1, -1, s*C_dict['M3']['Cdb'])
        add_y(N_3, -1, s*C_dict['M4']['Cgs']); add_y(N_3, N_2, s*C_dict['M4']['Cgd'])#; add_y(N_2, -1, s*C_dict['M4']['Cdb'])
        
        # M5/M6 (關鍵前饋路徑)
        add_y(N_1, -1, s*C_dict['M5']['Cgs']); add_y(N_1, N_P2, s*C_dict['M5']['Cgd']); add_y(N_P2, -1, s*C_dict['M5']['Cdb'])
        add_y(N_2, -1, s*C_dict['M6']['Cgs']); add_y(N_2, N_N2, s*C_dict['M6']['Cgd']); add_y(N_N2, -1, s*C_dict['M6']['Cdb'])
        
        # M7/M8
        add_y(N_4, -1, s*C_dict['M7']['Cgs']); add_y(N_4, N_P2, s*C_dict['M7']['Cgd']); add_y(N_P2, -1, s*C_dict['M7']['Cdb'])
        add_y(N_4, -1, s*C_dict['M8']['Cgs']); add_y(N_4, N_N2, s*C_dict['M8']['Cgd']); add_y(N_N2, -1, s*C_dict['M8']['Cdb'])
        
        # 尾節點對地總電容
        add_y(N_S, -1, s*(C_dict['M10']['Cdb'])) #  + C_WELL_S

        # 主動元件 gm, gds
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
    
    # ================= 虛擬失真電流計算 =================
    def calc_iNL2(K, vg, vd, mos_type):
        i_nl2 = K['K20']*(vg**2) + K['K11']*(vg*vd) + K['K02']*(vd**2)
        return -i_nl2 if mos_type == 'pch' else i_nl2

    def calc_iNL3(K, vg1, vd1, vg2, vd2, mos_type):
        i_pure3 = K['K30']*(vg1**3) + K['K21']*(vg1**2)*vd1 + K['K12']*vg1*(vd1**2) + K['K03']*(vd1**3)
        i_mix = 2*K['K20']*vg1*vg2 + K['K11']*(vg1*vd2 + vg2*vd1) + 2*K['K02']*vd1*vd2
        if mos_type == 'pch': return i_pure3 - i_mix
        else: return i_pure3 + i_mix

    def inject_inl(I_vec, d, src, inl_val):
        if d >= 0: I_vec[d] -= inl_val
        if src >= 0: I_vec[src] += inl_val

    # 3. 求解核心 (與之前邏輯相同)
    w = 2 * np.pi * FIN
    Y1 = build_Y_matrix(w)
    I_in1 = np.zeros(9, dtype=complex)
    I_in1[N_P1] = (0.5 * V_AMP_IN) / passives['R_in']
    I_in1[N_N1] = (-0.5 * V_AMP_IN) / passives['R_in']
    V1 = np.linalg.solve(Y1, I_in1)

    def get_vg_vd(V):
        return {
            'M1': (V[N_P1]-V[N_S], V[N_1]-V[N_S]),
            'M2': (V[N_N1]-V[N_S], V[N_2]-V[N_S]),
            'M3': (V[N_3], V[N_1]),
            'M4': (V[N_3], V[N_2]),
            'M5': (V[N_1], V[N_P2]),
            'M6': (V[N_2], V[N_N2]),
            'M7': (V[N_4], V[N_P2]),
            'M8': (V[N_4], V[N_N2]),
            'M10': (0, V[N_S])
        }
    vg1_dict = {k: v[0] for k, v in get_vg_vd(V1).items()}
    vd1_dict = {k: v[1] for k, v in get_vg_vd(V1).items()}

    I_in2 = np.zeros(9, dtype=complex)
    node_map = {'M1':(N_1, N_S), 'M2':(N_2, N_S),
                'M3':(N_1, -1), 'M4':(N_2, -1), 
                'M5':(N_P2, -1), 'M6':(N_N2, -1),
                'M7':(N_P2, -1), 'M8':(N_N2, -1),
                'M10':(N_S, -1)
    }
    
    for inst in node_map:
        inl2 = calc_iNL2(K_dict[inst], vg1_dict[inst], vd1_dict[inst], op_config[inst]['type_inl'])
        inject_inl(I_in2, node_map[inst][0], node_map[inst][1], inl2)
    
    Y2 = build_Y_matrix(2 * w)
    V2 = np.linalg.solve(Y2, I_in2)
    
    vg2_dict = {k: v[0] for k, v in get_vg_vd(V2).items()}
    vd2_dict = {k: v[1] for k, v in get_vg_vd(V2).items()}

    I_in3 = np.zeros(9, dtype=complex)
    for inst in node_map:
        inl3 = calc_iNL3(K_dict[inst], vg1_dict[inst], vd1_dict[inst], vg2_dict[inst], vd2_dict[inst], op_config[inst]['type_inl'])
        inject_inl(I_in3, node_map[inst][0], node_map[inst][1], inl3)
        
    Y3 = build_Y_matrix(3 * w)
    V3 = np.linalg.solve(Y3, I_in3)

    # 4. 指標結算
    v1_out = V1[N_P2] - V1[N_N2]
    v2_out = V2[N_P2] - V2[N_N2]
    v3_out = V3[N_P2] - V3[N_N2]

    a1_closed = v1_out / V_AMP_IN
    hd2_lin = 0.5 * np.abs(v2_out / v1_out)
    hd3_lin = 0.25 * np.abs(v3_out / v1_out)
    
    # print(f"\n=== Razavi 全差分 2-Stage OTA 閉環靜態失真 ===")
    # print(f"閉環差分增益 (A1): {np.abs(a1_closed):.5f} V/V")

    hd2_db = 20 * np.log10(hd2_lin) if hd2_lin > 1e-15 else -300
    hd3_db = 20 * np.log10(hd3_lin)

    # ================= 新增：OIP3 / IIP3 計算 =================
    # 1. 理論雙音互調失真 (IMD3 比 HD3 大 3 倍，即 9.54 dB)
    imd3_lin = 3.0 * hd3_lin
    imd3_db = 20 * np.log10(imd3_lin)
    
    # 2. 計算輸出基頻功率 (將差分輸出峰值振幅轉為 RMS，並算入 50 歐姆系統功率)
    R_SYS = 50.0
    vout_amp = np.abs(v1_out)          # 輸出峰值電壓振幅
    vout_rms = vout_amp / np.sqrt(2)   # 轉換為有效值
    
    # 計算基頻輸出功率 (dBm) 與 電壓位準 (dBV)
    pout_dbm = 10 * np.log10((vout_rms**2) / R_SYS / 1e-3)
    pout_dbv = 20 * np.log10(vout_rms)
    
    # 3. 截點公式：OIP3 = P_out + |IMD3| / 2
    # 注意：imd3_db 是負數 (例如 -80dBc)，所以我們用減法來加上它的絕對值
    oip3_dbm = pout_dbm - (imd3_db / 2.0)
    # oip3_dbv = pout_dbv - (imd3_db / 2.0)
    
    # 計算輸入三階截點 (IIP3) = OIP3 - 系統增益 (dB)
    gain_db = 20 * np.log10(np.abs(a1_closed))
    iip3_dbm = oip3_dbm - gain_db

    # ==========================================================
    return {
        'Gain_V_V': np.abs(a1_closed),
        'HD2_dBc': hd2_db,
        'HD3_dBc': hd3_db,
        'SFDR_dBc': -hd3_db,
        'IMD3_dBc': imd3_db,    # 附帶輸出雙音 IMD3 預測值
        'OIP3_dBm': oip3_dbm,   # 射頻工程師最愛的 OIP3
        'IIP3_dBm': iip3_dbm
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
    metrics = evaluate_volterra_sfdr(op_config, V_AMP_IN=0.316, FIN=53/512*100e6)
    
    print(f"HD2        : {metrics['HD2_dBc']:.2f} dBc (完美對稱抵消)")
    print(f"HD3        : {metrics['HD3_dBc']:.2f} dBc")
    print(f"SFDR       : {-metrics['HD3_dBc']:.2f} dBc")