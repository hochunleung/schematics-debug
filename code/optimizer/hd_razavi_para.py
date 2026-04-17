import numpy as np
import pandas as pd
import re
from lut_engine import MosData
from scipy.interpolate import RegularGridInterpolator

# ================= 核心求解流程 =================
def evaluate_volterra_sfdr(op_config, V_AMP_IN=0.632, FIN=10.35e6, F_SPACING=1e6):
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

    def calc_iNL2_mix(K, vgA, vdA, vgB, vdB, coeff, mos_type):
        """處理單音與雙音混頻的通用二階電流"""
        i_nl2 = coeff * (K['K20']*vgA*vgB + 0.5*K['K11']*(vgA*vdB + vgB*vdA) + K['K02']*vdA*vdB)
        return -i_nl2 if mos_type == 'pch' else i_nl2
    
    def calc_iNL3_imd3(K, va, vda, vb_c, vdb_c, v2a, vd2a, vdw, vddw, mos_type):
        """雙音 IMD3 專用三階電流 (精準捕捉低頻與高頻記憶效應)"""
        # 純三階互調 (展開後自帶係數 3)
        i_pure3 = 3*K['K30']*(va**2)*vb_c + K['K21']*(va**2*vdb_c + 2*va*vb_c*vda) + \
                  K['K12']*(vb_c*vda**2 + 2*va*vda*vdb_c) + 3*K['K03']*(vda**2)*vdb_c
        
        # 高頻混頻 (2*w1 與 -w2) & 低頻混頻 (w1 - w2 與 w1)
        i_mix_high = 2*K['K20']*v2a*vb_c + K['K11']*(v2a*vdb_c + vb_c*vd2a) + 2*K['K02']*vd2a*vdb_c
        i_mix_low = 2*K['K20']*vdw*va + K['K11']*(vdw*vda + va*vddw) + 2*K['K02']*vddw*vda
        
        i_mix = i_mix_high + i_mix_low
        return i_pure3 - i_mix if mos_type == 'pch' else i_pure3 + i_mix
    
    def inject_inl(I_vec, d, src, inl_val):
        if d >= 0: I_vec[d] -= inl_val
        if src >= 0: I_vec[src] += inl_val

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
    
    node_map = {'M1':(N_1, N_S), 'M2':(N_2, N_S),
                'M3':(N_1, -1), 'M4':(N_2, -1), 
                'M5':(N_P2, -1), 'M6':(N_N2, -1),
                'M7':(N_P2, -1), 'M8':(N_N2, -1),
                'M10':(N_S, -1)
    }

    # 3. 求解核心 (與之前邏輯相同)
    # 定義所有需要的頻率點
    w = 2 * np.pi * FIN
    w1 = 2 * np.pi * FIN
    w2 = 2 * np.pi * (FIN + F_SPACING)
    w_2a = 2 * w1          # 二次諧波 (供 HD2, HD3 混頻用)
    w_dw = w1 - w2         # 低頻包絡 (供 IMD3 混頻用)
    w_3a = 3 * w1          # 單音 HD3 輸出頻率
    w_imd3 = 2 * w1 - w2   # 雙音 IMD3 輸出頻率

    # [一階求解]
    Y1a = build_Y_matrix(w1)
    Y1b = build_Y_matrix(w2)
    I_in1 = np.zeros(9, dtype=complex)
    I_in1[N_P1] = (0.5 * V_AMP_IN) / passives['R_in']
    I_in1[N_N1] = (-0.5 * V_AMP_IN) / passives['R_in']

    V1a = np.linalg.solve(Y1a, I_in1)
    V1b = np.linalg.solve(Y1b, I_in1)
    #vg1a, vd1a = get_vg_vd(V1a)
    vg1a = {k: v[0] for k, v in get_vg_vd(V1a).items()}
    vd1a = {k: v[1] for k, v in get_vg_vd(V1a).items()}
    #vg1b, vd1b = get_vg_vd(V1b)
    vg1b = {k: v[0] for k, v in get_vg_vd(V1b).items()}
    vd1b = {k: v[1] for k, v in get_vg_vd(V1b).items()}

    # [二階求解]
    Y2a = build_Y_matrix(w_2a)
    Ydw = build_Y_matrix(w_dw)
    I_in2_2a = np.zeros(9, dtype=complex)
    I_in2_dw = np.zeros(9, dtype=complex)

    for inst in node_map:
        K, mt = K_dict[inst], op_config[inst]['type_inl']
        # 1. 自身平方 (係數 1.0)
        inl2_2a = calc_iNL2_mix(K, vg1a[inst], vd1a[inst], vg1a[inst], vd1a[inst], 1.0, mt)
        # 2. 差頻混頻 (w1 與 -w2, 係數 2.0, w2 必須取共軛)
        inl2_dw = calc_iNL2_mix(K, vg1a[inst], vd1a[inst], np.conj(vg1b[inst]), np.conj(vd1b[inst]), 2.0, mt)
        
        inject_inl(I_in2_2a, node_map[inst][0], node_map[inst][1], inl2_2a)
        inject_inl(I_in2_dw, node_map[inst][0], node_map[inst][1], inl2_dw)
    
    V2a = np.linalg.solve(Y2a, I_in2_2a)
    Vdw = np.linalg.solve(Ydw, I_in2_dw)

    #vg2a, vd2a = get_vg_vd(V2a)
    vg2a = {k: v[0] for k, v in get_vg_vd(V2a).items()}
    vd2a = {k: v[1] for k, v in get_vg_vd(V2a).items()}
    #vgdw, vddw = get_vg_vd(Vdw)
    vgdw = {k: v[0] for k, v in get_vg_vd(Vdw).items()}
    vddw = {k: v[1] for k, v in get_vg_vd(Vdw).items()}

    # [三階求解]
    Y3a = build_Y_matrix(w_3a)
    Y3imd3 = build_Y_matrix(w_imd3)
    I_in3_3a = np.zeros(9, dtype=complex)
    I_in3_imd3 = np.zeros(9, dtype=complex)
    
    for inst in node_map:
        K, mt = K_dict[inst], op_config[inst]['type_inl']
        
        # 1. 單音 HD3 激勵
        inl3_3a = calc_iNL3(K, vg1a[inst], vd1a[inst], vg2a[inst], vd2a[inst], mt)
        inject_inl(I_in3_3a, node_map[inst][0], node_map[inst][1], inl3_3a)
        
        # 2. 雙音 IMD3 激勵 (精準記憶效應注入)
        inl3_imd3 = calc_iNL3_imd3(K, vg1a[inst], vd1a[inst], np.conj(vg1b[inst]), np.conj(vd1b[inst]), 
                                   vg2a[inst], vd2a[inst], vgdw[inst], vddw[inst], mt)
        inject_inl(I_in3_imd3, node_map[inst][0], node_map[inst][1], inl3_imd3)
        
    V3a = np.linalg.solve(Y3a, I_in3_3a)
    V3_imd3 = np.linalg.solve(Y3imd3, I_in3_imd3)

    # 4. 指標結算
    v1_out = V1a[N_P2] - V1a[N_N2]
    v2_out = V2a[N_P2] - V2a[N_N2]
    v3_out_hd3 = V3a[N_P2] - V3a[N_N2]
    v3_out_imd3 = V3_imd3[N_P2] - V3_imd3[N_N2]

    # 線性相量振幅轉化 (數學上都需要乘上 0.25 來對應時域 Cosine)
    a1_closed = v1_out / V_AMP_IN
    hd2_lin = 0.25 * np.abs(v2_out / v1_out)  # 雙輸入 0.5 * 0.5
    hd3_lin = 0.25 * np.abs(v3_out_hd3 / v1_out)
    imd3_lin = 0.25 * np.abs(v3_out_imd3 / v1_out)
    
    hd2_db = 20 * np.log10(hd2_lin) if hd2_lin > 1e-15 else -300
    hd3_db = 20 * np.log10(hd3_lin)
    imd3_db = 20 * np.log10(imd3_lin)

    # 射頻 OIP3 結算 (基於 50 歐姆系統)
    R_SYS = 50.0
    vout_rms = np.abs(v1_out) / np.sqrt(2)
    pout_dbm = 10 * np.log10((vout_rms**2) / R_SYS / 1e-3)
    oip3_dbm = pout_dbm - (imd3_db / 2.0)
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