import numpy as np
import pandas as pd
import re
from lut_engine import MosData
from scipy.interpolate import RegularGridInterpolator

# ================= 核心求解流程 =================
def evaluate_volterra_sfdr(op_config, V_AMP_IN=0.632, FIN=10.35e6, F_SPACING=1e6):
    # ================= 節點索引定義 (擴充至 13 維 MNA) =================
    N_P1, N_N1, N_S = 0, 1, 2
    N_1, N_2, N_3 = 3, 4, 5
    N_P2, N_N2, N_4 = 6, 7, 8
    N_INP, N_INN = 9, 10       # 新增：真正的外部輸入節點
    I_VINP, I_VINN = 11, 12    # 新增：MNA 引入的理想電壓源電流
    # ====================================================================

    # 1. 動態提取 K 字典和寄生電容字典
    K_dict = {}
    C_dict = {}
    for inst, cfg in op_config['mosfets'].items():
        m = cfg['model']
        L, W, VGS, VDS = cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
        
        K_dict[inst] = {k: m.lookup_by_vgs(k+'_W', L, VGS, VDS) * W for k in ['K10', 'K01', 'K20', 'K02', 'K11', 'K30', 'K03', 'K21', 'K12']}
        C_dict[inst] = {
            'Cgs': m.lookup_by_vgs('CGS_W', L, VGS, VDS) * W,
            'Cgd': m.lookup_by_vgs('CGD_W', L, VGS, VDS) * W,
            'Cdb': np.abs(m.lookup_by_vgs('CDB_W', L, VGS, VDS) * W) 
        }

    passives = op_config['passives']

    # 2. 構建動態 Y 矩陣 (13x13 MNA)
    def build_Y_matrix(w):
        Y = np.zeros((13, 13), dtype=complex)
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
        
        # 負回饋電阻 (R34, R35)
        add_y(N_P1, N_N2, 1/R_in); add_y(N_N1, N_P2, 1/R_in)
        
        # 關鍵修改：輸入電阻 R37, R36 連接到真正的 N_INP 和 N_INN 節點！
        add_y(N_P1, N_INP, 1/R_in); add_y(N_N1, N_INN, 1/R_in)
        
        # === MNA 擴充矩陣核心：加入理想電壓源的 KVL 邊界條件 ===
        Y[N_INP, I_VINP] = 1.0; Y[I_VINP, N_INP] = 1.0
        Y[N_INN, I_VINN] = 1.0; Y[I_VINN, N_INN] = 1.0
        # =======================================================
        
        add_y(N_P2, N_4, 1/R_cmfb); add_y(N_N2, N_4, 1/R_cmfb)
        add_y(N_1, N_3, 1/R_cmfb); add_y(N_2, N_3, 1/R_cmfb)
        
        y_miller = 1.0 / (passives['Rz'] + 1.0/(s * passives['Cc']))
        add_y(N_P2, N_1, y_miller); add_y(N_N2, N_2, y_miller)

        # 動態電容注入
        add_y(N_P1, N_S, s*C_dict['M1']['Cgs']); add_y(N_P1, N_1, s*C_dict['M1']['Cgd'])
        add_y(N_N1, N_S, s*C_dict['M2']['Cgs']); add_y(N_N1, N_2, s*C_dict['M2']['Cgd'])
        add_y(N_3, -1, s*C_dict['M3']['Cgs']); add_y(N_3, N_1, s*C_dict['M3']['Cgd'])
        add_y(N_3, -1, s*C_dict['M4']['Cgs']); add_y(N_3, N_2, s*C_dict['M4']['Cgd'])
        add_y(N_1, -1, s*C_dict['M5']['Cgs']); add_y(N_1, N_P2, s*C_dict['M5']['Cgd']); add_y(N_P2, -1, s*C_dict['M5']['Cdb'])
        add_y(N_2, -1, s*C_dict['M6']['Cgs']); add_y(N_2, N_N2, s*C_dict['M6']['Cgd']); add_y(N_N2, -1, s*C_dict['M6']['Cdb'])
        add_y(N_4, -1, s*C_dict['M7']['Cgs']); add_y(N_4, N_P2, s*C_dict['M7']['Cgd']); add_y(N_P2, -1, s*C_dict['M7']['Cdb'])
        add_y(N_4, -1, s*C_dict['M8']['Cgs']); add_y(N_4, N_N2, s*C_dict['M8']['Cgd']); add_y(N_N2, -1, s*C_dict['M8']['Cdb'])
        
        # 尾節點對地總電容
        add_y(N_S, -1, s*(C_dict['M10']['Cdb'])) 

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
        i_nl2 = coeff * (K['K20']*vgA*vgB + 0.5*K['K11']*(vgA*vdB + vgB*vdA) + K['K02']*vdA*vdB)
        return -i_nl2 if mos_type == 'pch' else i_nl2
    
    def calc_iNL3_imd3(K, va, vda, vb_c, vdb_c, v2a, vd2a, vdw, vddw, mos_type):
        i_pure3 = 3*K['K30']*(va**2)*vb_c + K['K21']*(va**2*vdb_c + 2*va*vb_c*vda) + \
                  K['K12']*(vb_c*vda**2 + 2*va*vda*vdb_c) + 3*K['K03']*(vda**2)*vdb_c
        
        i_mix_high = 2*K['K20']*v2a*vb_c + K['K11']*(v2a*vdb_c + vb_c*vd2a) + 2*K['K02']*vd2a*vdb_c
        i_mix_low = 2*K['K20']*vdw*va + K['K11']*(vdw*vda + va*vddw) + 2*K['K02']*vddw*vda
        
        i_mix = i_mix_high + i_mix_low
        return i_pure3 - i_mix if mos_type == 'pch' else i_pure3 + i_mix
    
    def inject_inl(I_vec, d, src, inl_val):
        if d >= 0: I_vec[d] -= inl_val
        if src >= 0: I_vec[src] += inl_val

    def get_vg_vd(V):
        return {
            'M1': (V[N_P1]-V[N_S], V[N_1]-V[N_S]), 'M2': (V[N_N1]-V[N_S], V[N_2]-V[N_S]),
            'M3': (V[N_3], V[N_1]), 'M4': (V[N_3], V[N_2]),
            'M5': (V[N_1], V[N_P2]), 'M6': (V[N_2], V[N_N2]),
            'M7': (V[N_4], V[N_P2]), 'M8': (V[N_4], V[N_N2]),
            'M10': (0, V[N_S])
        }
    
    node_map = {'M1':(N_1, N_S), 'M2':(N_2, N_S),
                'M3':(N_1, -1), 'M4':(N_2, -1), 
                'M5':(N_P2, -1), 'M6':(N_N2, -1),
                'M7':(N_P2, -1), 'M8':(N_N2, -1),
                'M10':(N_S, -1)}

    # 3. 求解核心
    w1 = 2 * np.pi * FIN
    w2 = 2 * np.pi * (FIN + F_SPACING)
    w_2a = 2 * w1          
    w_dw = w1 - w2         
    w_3a = 3 * w1          
    w_imd3 = 2 * w1 - w2   

    # [一階求解]
    Y1a = build_Y_matrix(w1)
    Y1b = build_Y_matrix(w2)
    
    # 建立 13 維激勵向量，並在 MNA 的電壓變數列直接注入理想電壓！
    I_in1 = np.zeros(13, dtype=complex)
    I_in1[I_VINP] = 0.5 * V_AMP_IN
    I_in1[I_VINN] = -0.5 * V_AMP_IN

    V1a = np.linalg.solve(Y1a, I_in1)
    V1b = np.linalg.solve(Y1b, I_in1)
    
    vg1a = {k: v[0] for k, v in get_vg_vd(V1a).items()}
    vd1a = {k: v[1] for k, v in get_vg_vd(V1a).items()}
    vg1b = {k: v[0] for k, v in get_vg_vd(V1b).items()}
    vd1b = {k: v[1] for k, v in get_vg_vd(V1b).items()}

    # [二階求解]
    Y2a = build_Y_matrix(w_2a)
    Ydw = build_Y_matrix(w_dw)
    I_in2_2a = np.zeros(13, dtype=complex)
    I_in2_dw = np.zeros(13, dtype=complex)

    for inst in node_map:
        K, mt = K_dict[inst], op_config['mosfets'][inst]['type_inl']
        inl2_2a = calc_iNL2_mix(K, vg1a[inst], vd1a[inst], vg1a[inst], vd1a[inst], 1.0, mt)
        inl2_dw = calc_iNL2_mix(K, vg1a[inst], vd1a[inst], np.conj(vg1b[inst]), np.conj(vd1b[inst]), 2.0, mt)
        
        inject_inl(I_in2_2a, node_map[inst][0], node_map[inst][1], inl2_2a)
        inject_inl(I_in2_dw, node_map[inst][0], node_map[inst][1], inl2_dw)
    
    V2a = np.linalg.solve(Y2a, I_in2_2a)
    Vdw = np.linalg.solve(Ydw, I_in2_dw)

    vg2a = {k: v[0] for k, v in get_vg_vd(V2a).items()}
    vd2a = {k: v[1] for k, v in get_vg_vd(V2a).items()}
    vgdw = {k: v[0] for k, v in get_vg_vd(Vdw).items()}
    vddw = {k: v[1] for k, v in get_vg_vd(Vdw).items()}

    # [三階求解]
    Y3a = build_Y_matrix(w_3a)
    Y3imd3 = build_Y_matrix(w_imd3)
    I_in3_3a = np.zeros(13, dtype=complex)
    I_in3_imd3 = np.zeros(13, dtype=complex)
    
    for inst in node_map:
        K, mt = K_dict[inst], op_config['mosfets'][inst]['type_inl']
        
        inl3_3a = calc_iNL3(K, vg1a[inst], vd1a[inst], vg2a[inst], vd2a[inst], mt)
        inject_inl(I_in3_3a, node_map[inst][0], node_map[inst][1], inl3_3a)
        
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

    a1_closed = v1_out / V_AMP_IN
    hd2_lin = 0.25 * np.abs(v2_out / v1_out)  
    hd3_lin = 0.25 * np.abs(v3_out_hd3 / v1_out)
    imd3_lin = 0.25 * np.abs(v3_out_imd3 / v1_out)
    
    hd2_db = 20 * np.log10(hd2_lin) if hd2_lin > 1e-15 else -300
    hd3_db = 20 * np.log10(hd3_lin)
    imd3_db = 20 * np.log10(imd3_lin)

    R_SYS = 50.0
    vout_rms = np.abs(v1_out) / np.sqrt(2)
    pout_dbm = 10 * np.log10((vout_rms**2) / R_SYS / 1e-3)
    oip3_dbm = pout_dbm - (imd3_db / 2.0)
    
    gain_db = 20 * np.log10(np.abs(a1_closed))
    iip3_dbm = oip3_dbm - gain_db

    return {
        'Gain_V_V': np.abs(a1_closed),
        'HD2_dBc': hd2_db,
        'HD3_dBc': hd3_db,
        'SFDR_dBc': -hd3_db,
        'IMD3_dBc': imd3_db,
        'OIP3_dBm': oip3_dbm,
        'IIP3_dBm': iip3_dbm
    }


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
        },
        # 被動元件一併傳遞
        'passives': {'R_f': 4000, 'R_in': 4000, 'Rz': 204.6, 'Cc': 402e-15, 'R_cmfb': 40e3}
    }
    metrics = evaluate_volterra_sfdr(op_config, V_AMP_IN=0.316, FIN=53/512*100e6)
    
    print(f"HD2        : {metrics['HD2_dBc']:.2f} dBc (完美對稱抵消)")
    print(f"HD3        : {metrics['HD3_dBc']:.2f} dBc")
    print(f"SFDR       : {-metrics['HD3_dBc']:.2f} dBc")