import numpy as np
import pandas as pd
import re
from lut_engine import MosData
from scipy.interpolate import RegularGridInterpolator
from pprint import pprint

def hd_analysis(op_config, V_AMP_IN=0.632, FIN=10.35e6, F_SPACING=1e6):
    """
    通用化的 Volterra 分析器，專為任意電路拓撲設計。
    :param op_config: 包含電路拓撲、操作點和元件參數的字典。
                      預期包含 'all_nodes', 'mosfets', 'R', 'C', 'series_RC', 'passives', 'hd_probes'。
    :param V_AMP_IN: 輸入信號的峰值振幅。
    :param FIN: 輸入信號的中心頻率。
    :param F_SPACING: 雙音測試時的頻率間隔。
    :return: dict 包含 Gain_V_V, HD2_dBc, HD3_dBc, SFDR_dBc, IMD3_dBc, OIP3_dBm, IIP3_dBm
    """
    # =================================================================
    # 1. 矩陣索引映射 (Mapping)
    # =================================================================
    GND = -1 # AC 地的約定索引

    all_nodes = op_config['all_nodes']
    # 過濾掉固定電壓節點，這些節點不作為 MNA 矩陣的變量
    # '0' (GND) 通過 GND = -1 處理，因此不應包含在 dynamic_nodes 中
    fixed_voltages_to_exclude = ['VDD', 'VCM', '0'] 
    
    # 創建將包含在 MNA 矩陣中的動態節點列表
    dynamic_nodes = sorted(list(set(all_nodes) - set(fixed_voltages_to_exclude)))
    
    NUM_VARS = len(dynamic_nodes)
    node_to_idx = {n: i for i, n in enumerate(dynamic_nodes)}
    
    # 1. 動態提取 K 字典和寄生電容字典
    K_dict = {}
    C_dict = {}
    for inst, cfg in op_config['mosfets'].items():
        m = cfg['model']
        L, W, VGS, VDS = cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
        
        # 提取 K 參數 (密度 * 實際寬度)
        K_dict[inst] = {k: m.lookup_by_vgs(k+'_W', L, VGS, VDS) * W for k in ['K10', 'K01', 'K20', 'K02', 'K11', 'K30', 'K03', 'K21', 'K12']}
        
        # 提取寄生電容 (密度 * 實際寬度，並確保 Cdb 為正絕對值)
        C_dict[inst] = {
            'Cgs': m.lookup_by_vgs('CGS_W', L, VGS, VDS) * W,
            'Cgd': m.lookup_by_vgs('CGD_W', L, VGS, VDS) * W,
            'Cdb': np.abs(m.lookup_by_vgs('CDB_W', L, VGS, VDS) * W) 
        }

    # 2. 構建動態 Y 矩陣
    def build_Y_matrix(w):
        Y = np.zeros((NUM_VARS, NUM_VARS), dtype=complex)
        s = 1j * w
        
        def add_y(name1, name2, y_val):
            n1 = node_to_idx[name1] if name1 in node_to_idx else GND
            n2 = node_to_idx[name2] if name2 in node_to_idx else GND
            if n1 != GND: Y[n1, n1] += y_val
            if n2 != GND: Y[n2, n2] += y_val
            if n1 != GND and n2 != GND:
                Y[n1, n2] -= y_val; Y[n2, n1] -= y_val

        def add_gm(d_name, g_name, src_name, gm):
            d = node_to_idx[d_name] if d_name in node_to_idx else GND
            g = node_to_idx[g_name] if g_name in node_to_idx else GND
            src = node_to_idx[src_name] if src_name in node_to_idx else GND
            """添加 VCCS (跨導)"""
            if d != GND:
                if g != GND:   Y[d, g] += gm
                if src != GND: Y[d, src] -= gm
            if src != GND:
                if g != GND:   Y[src, g] -= gm
                if src != GND: Y[src, src] += gm

        # 被動元件
        for res_name, res_cfg in op_config.get('R', {}).items():
            add_y(res_cfg['p'], res_cfg['n'], 1/res_cfg['val'])
        for cap_name, cap_cfg in op_config.get('C', {}).items():
            add_y(cap_cfg['p'], cap_cfg['n'], s * cap_cfg['val'])
        
        # 串聯 R-C 元件 (例如米勒補償)
        # for rc_name, rc_cfg in op_config.get('series_RC', {}).items():
        #     y_series_rc = (s * rc_cfg['C']) / (1 + s * rc_cfg['R'] * rc_cfg['C'])
        #     add_y(rc_cfg['p'], rc_cfg['n'], y_series_rc)

        # 主動元件 gm, gds, 寄生電容
        for mos_name, mos_cfg in op_config['mosfets'].items():
            d_node = mos_cfg['d']
            g_node = mos_cfg['g']
            s_node = mos_cfg['s']
            b_node = mos_cfg['b'] # 體節點，通常連接到源極或地

            # gm 和 gds
            add_gm(d_node, g_node, s_node, K_dict[mos_name]['K10'])
            add_y(d_node, s_node, K_dict[mos_name]['K01']) # gds 在 D 和 S 之間

            # 寄生電容
            add_y(g_node, s_node, s * C_dict[mos_name]['Cgs'])
            add_y(d_node, g_node, s * C_dict[mos_name]['Cgd'])
            add_y(d_node, b_node, s * C_dict[mos_name]['Cdb']) # Cdb 在 D 和 B 之間

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
    
    def inject_inl(I_vec, d_name, s_name, inl_val):
        d = node_to_idx[d_name] if d_name in node_to_idx else GND
        s = node_to_idx[s_name] if s_name in node_to_idx else GND
        if d != GND: I_vec[d] -= inl_val
        if s != GND: I_vec[s] += inl_val

    def get_vg_vd_from_dict(V_solution_dict):
        vg_vd_dict = {}
        for inst, cfg in op_config['mosfets'].items():
            g_node_name = cfg['g']
            s_node_name = cfg['s']
            d_node_name = cfg['d']

            vg = V_solution_dict.get(g_node_name, 0)
            vs = V_solution_dict.get(s_node_name, 0)
            vd = V_solution_dict.get(d_node_name, 0)

            vg_eff = vg - vs # Vgs
            vd_eff = vd - vs # Vds
            vg_vd_dict[inst] = (vg_eff, vd_eff)
        return vg_vd_dict
    
    # 3. 求解核心
    w1 = 2 * np.pi * FIN
    w2 = 2 * np.pi * (FIN + F_SPACING)
    w_2a = 2 * w1          # 二次諧波 (供 HD2, HD3 混頻用)
    w_dw = w1 - w2         # 低頻包絡 (供 IMD3 混頻用)
    w_3a = 3 * w1          # 單音 HD3 輸出頻率
    w_imd3 = 2 * w1 - w2   # 雙音 IMD3 輸出頻率

    # [一階求解]
    Y1a = build_Y_matrix(w1)
    Y1b = build_Y_matrix(w2)
    I_in1 = np.zeros(NUM_VARS, dtype=complex)
    
    input_p_node = op_config['hd_probes']['input']['p']
    input_n_node = op_config['hd_probes']['input']['n']
    
    R_in_val = op_config['passives']['R_in'] # 假設 R_in 始終在 passives 中用於輸入電流計算

    I_in1[node_to_idx[input_p_node]] = (0.5 * V_AMP_IN) / R_in_val
    I_in1[node_to_idx[input_n_node]] = (-0.5 * V_AMP_IN) / R_in_val

    V1a_raw = np.linalg.solve(Y1a, I_in1)
    V1b_raw = np.linalg.solve(Y1b, I_in1)

    # 將原始解轉換為節點名稱到電壓的字典
    V1a = {node: V1a_raw[idx] for node, idx in node_to_idx.items()}
    V1b = {node: V1b_raw[idx] for node, idx in node_to_idx.items()}
    
    vg1a_tuple = get_vg_vd_from_dict(V1a)
    vg1a = {k: v[0] for k, v in vg1a_tuple.items()}
    vd1a = {k: v[1] for k, v in vg1a_tuple.items()}
    
    vg1b_tuple = get_vg_vd_from_dict(V1b)
    vg1b = {k: v[0] for k, v in vg1b_tuple.items()}
    vd1b = {k: v[1] for k, v in vg1b_tuple.items()}

    # [二階求解]
    Y2a = build_Y_matrix(w_2a)
    Ydw = build_Y_matrix(w_dw)
    I_in2_2a = np.zeros(NUM_VARS, dtype=complex)
    I_in2_dw = np.zeros(NUM_VARS, dtype=complex)
    
    for inst, mos_cfg in op_config['mosfets'].items():
        K, mt = K_dict[inst], mos_cfg['type_inl'] # 假設 'type_inl' 在 mos_cfg 中
        d_node_name = mos_cfg['d']
        s_node_name = mos_cfg['s']

        # 1. 自身平方 (係數 1.0)
        inl2_2a = calc_iNL2_mix(K, vg1a[inst], vd1a[inst], vg1a[inst], vd1a[inst], 1.0, mt)
        # 2. 差頻混頻 (w1 與 -w2, 係數 2.0, w2 必須取共軛)
        inl2_dw = calc_iNL2_mix(K, vg1a[inst], vd1a[inst], np.conj(vg1b[inst]), np.conj(vd1b[inst]), 2.0, mt)
        
        inject_inl(I_in2_2a, d_node_name, s_node_name, inl2_2a)
        inject_inl(I_in2_dw, d_node_name, s_node_name, inl2_dw)
    
    V2a_raw = np.linalg.solve(Y2a, I_in2_2a)
    Vdw_raw = np.linalg.solve(Ydw, I_in2_dw)

    V2a = {node: V2a_raw[idx] for node, idx in node_to_idx.items()}
    Vdw = {node: Vdw_raw[idx] for node, idx in node_to_idx.items()}

    vg2a_tuple = get_vg_vd_from_dict(V2a)
    vg2a = {k: v[0] for k, v in vg2a_tuple.items()}
    vd2a = {k: v[1] for k, v in vg2a_tuple.items()}
    
    vgdw_tuple = get_vg_vd_from_dict(Vdw)
    vgdw = {k: v[0] for k, v in vgdw_tuple.items()}
    vddw = {k: v[1] for k, v in vgdw_tuple.items()}

    # [三階求解]
    Y3a = build_Y_matrix(w_3a)
    Y3imd3 = build_Y_matrix(w_imd3)
    I_in3_3a = np.zeros(NUM_VARS, dtype=complex)
    I_in3_imd3 = np.zeros(NUM_VARS, dtype=complex)
    
    for inst, mos_cfg in op_config['mosfets'].items():
        K, mt = K_dict[inst], mos_cfg['type_inl']
        d_node_name = mos_cfg['d']
        s_node_name = mos_cfg['s']
        
        # 1. 單音 HD3 激勵
        inl3_3a = calc_iNL3(K, vg1a[inst], vd1a[inst], vg2a[inst], vd2a[inst], mt)
        inject_inl(I_in3_3a, d_node_name, s_node_name, inl3_3a)
        
        # 2. 雙音 IMD3 激勵 (精準記憶效應注入)
        inl3_imd3 = calc_iNL3_imd3(K, vg1a[inst], vd1a[inst], np.conj(vg1b[inst]), np.conj(vd1b[inst]), 
                                   vg2a[inst], vd2a[inst], vgdw[inst], vddw[inst], mt)
        inject_inl(I_in3_imd3, d_node_name, s_node_name, inl3_imd3)
        
    V3a_raw = np.linalg.solve(Y3a, I_in3_3a)
    V3_imd3_raw = np.linalg.solve(Y3imd3, I_in3_imd3)

    V3a = {node: V3a_raw[idx] for node, idx in node_to_idx.items()}
    V3_imd3 = {node: V3_imd3_raw[idx] for node, idx in node_to_idx.items()}

    # 4. 指標結算
    output_p_node = op_config['hd_probes']['output']['p']
    output_n_node = op_config['hd_probes']['output']['n']

    v1_out = V1a.get(output_p_node, 0) - V1a.get(output_n_node, 0)
    v2_out = V2a.get(output_p_node, 0) - V2a.get(output_n_node, 0)
    v3_out_hd3 = V3a.get(output_p_node, 0) - V3a.get(output_n_node, 0)
    v3_out_imd3 = V3_imd3.get(output_p_node, 0) - V3_imd3.get(output_n_node, 0)

    # 線性相量振幅轉化 (數學上都需要乘上 0.25 來對應時域 Cosine)
    a1_closed = v1_out / V_AMP_IN
    hd2_lin = 0.25 * np.abs(v2_out / v1_out)  # 雙輸入 0.5 * 0.5
    hd3_lin = 0.25 * np.abs(v3_out_hd3 / v1_out)
    imd3_lin = 0.25 * np.abs(v3_out_imd3 / v1_out)
    
    hd2_db = 20 * np.log10(hd2_lin) if hd2_lin > 1e-15 else -300
    hd3_db = 20 * np.log10(hd3_lin) if hd3_lin > 1e-15 else -300
    imd3_db = 20 * np.log10(imd3_lin) if imd3_lin > 1e-15 else -300

    # 射頻 OIP3 結算 (基於 50 歐姆系統)
    R_SYS = 50.0
    vout_rms = np.abs(v1_out) / np.sqrt(2)
    # 處理 vout_rms 過小導致的除零或對數為零問題
    if vout_rms > 1e-15:
        pout_dbm = 10 * np.log10((vout_rms**2) / R_SYS / 1e-3)
    else:
        pout_dbm = -300 # 非常小的功率
    
    oip3_dbm = pout_dbm - (imd3_db / 2.0)
    # 計算輸入三階截點 (IIP3) = OIP3 - 系統增益 (dB)
    gain_db = 20 * np.log10(np.abs(a1_closed)) if np.abs(a1_closed) > 1e-15 else -300
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

# ================= 測試代碼 (需要一個 op_config 範例) =================
if __name__ == "__main__":
    # 載入 LUT 資料庫 (注意 P 管開啟 is_pmos=True)
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
    models = {'nch': nch, 'nch_lvt': nch_lvt, 'pch': pch}

    # 模擬一個 op_config，基於 hd_razavi_para.py 中的例子
    test_op_config = {
        # 必須包含 mid1 和 mid2，否則矩陣會將其視為 GND
        'all_nodes': ['N_P1', 'N_N1', 'N_S', 'N_1', 'N_2', 'N_3', 'N_P2', 'N_N2', 'N_4', 'mid1', 'mid2', 'VDD', 'VCM', '0'],
        'mosfets': {
            'M1':  {'d': 'N_1', 'g': 'N_P1', 's': 'N_S', 'b': 'N_S', 'VGS': 0.2565, 'VDS': 0.4191, 'W': 96e-6, 'L': 100e-9, 'model': nch_lvt, 'type_inl': 'nch'},
            'M2':  {'d': 'N_2', 'g': 'N_N1', 's': 'N_S', 'b': 'N_S', 'VGS': 0.2565, 'VDS': 0.4191, 'W': 96e-6, 'L': 100e-9, 'model': nch_lvt, 'type_inl': 'nch'},
            'M3':  {'d': 'N_1', 'g': 'N_3', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.2873, 'W': 96e-6, 'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
            'M4':  {'d': 'N_2', 'g': 'N_3', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.2873, 'W': 96e-6, 'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
            'M5':  {'d': 'N_P2', 'g': 'N_1', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.4498, 'W': 192e-6,'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
            'M6':  {'d': 'N_N2', 'g': 'N_2', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.4498, 'W': 192e-6,'L': 100e-9, 'model': pch, 'type_inl': 'pch'},
            'M7':  {'d': 'N_P2', 'g': 'N_4', 's': '0', 'b': '0', 'VGS': 0.3301, 'VDS': 0.4502, 'W': 48e-6, 'L': 100e-9, 'model': nch, 'type_inl': 'nch'},
            'M8':  {'d': 'N_N2', 'g': 'N_4', 's': '0', 'b': '0', 'VGS': 0.3301, 'VDS': 0.4502, 'W': 48e-6, 'L': 100e-9, 'model': nch, 'type_inl': 'nch'},
            'M10': {'d': 'N_S', 'g': 'VCM', 's': '0', 'b': '0', 'VGS': 0.4585, 'VDS': 0.1936, 'W': 32e-6, 'L': 500e-9, 'model': nch, 'type_inl': 'nch'},
        },
        'R': {
            'R_f_p': {'p': 'N_P1', 'n': 'N_N2', 'val': 4000},
            'R_f_n': {'p': 'N_N1', 'n': 'N_P2', 'val': 4000},
            'R_in_p': {'p': 'N_P1', 'n': '0', 'val': 4000},
            'R_in_n': {'p': 'N_N1', 'n': '0', 'val': 4000},
            'R_cmfb_p2_4': {'p': 'N_P2', 'n': 'N_4', 'val': 40e3},
            'R_cmfb_n2_4': {'p': 'N_N2', 'n': 'N_4', 'val': 40e3},
            'R_cmfb_1_3': {'p': 'N_1', 'n': 'N_3', 'val': 40e3},
            'R_cmfb_2_3': {'p': 'N_2', 'n': 'N_3', 'val': 40e3},
            'Rz1': {'p': 'mid1', 'n': 'N_1', 'val': 204.6},
            'Rz2': {'p': 'mid2', 'n': 'N_2', 'val': 204.6},
        },
        'C': {
            # 在此範例中，沒有獨立的 Cc，它作為 series_RC 的一部分
            'Cc1': {'p': 'N_P2', 'n': 'mid1', 'val': 402e-15},
            'Cc2': {'p': 'N_N2', 'n': 'mid2', 'val': 402e-15},
        },
        # 'series_RC': {
        #     'miller_p': {'p': 'N_P2', 'n': 'N_1', 'R': 204.6, 'C': 402e-15},
        #     'miller_n': {'p': 'N_N2', 'n': 'N_2', 'R': 204.6, 'C': 402e-15},
        # },
        'passives': {'R_in': 4000}, # 此 R_in 用於輸入電流計算，不一定是 MNA 矩陣中的物理電阻。
        'hd_probes': {
            'input': {'p': 'N_P1', 'n': 'N_N1'},
            'output': {'p': 'N_P2', 'n': 'N_N2'}
        }
    }
    
    metrics = hd_analysis(test_op_config, V_AMP_IN=0.316, FIN=53/512*100e6)
    
    print(f"HD2        : {metrics['HD2_dBc']:.2f} dBc")
    print(f"HD3        : {metrics['HD3_dBc']:.2f} dBc")
    print(f"SFDR       : {-metrics['HD3_dBc']:.2f} dBc")
    print(f"IMD3       : {metrics['IMD3_dBc']:.2f} dBc")
    print(f"OIP3       : {metrics['OIP3_dBm']:.2f} dBm")
    print(f"IIP3       : {metrics['IIP3_dBm']:.2f} dBm")
