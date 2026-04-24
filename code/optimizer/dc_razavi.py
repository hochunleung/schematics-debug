import numpy as np
from lut_engine import MosData
from pprint import pprint

def solve_dc_fully_diff(param, models, vdd_val=0.9, vcm_val=0.45, W_unit=4e-6, L_unit=10e-9):
    nch = models['nch']
    nch_lvt = models['nch_lvt']
    pch = models['pch']
    dims = {
        'W1': param[0] * W_unit, 'L1': param[1] * L_unit,
        'W3': param[2] * W_unit, 'L3': param[3] * L_unit,
        'W5': param[4] * W_unit, 'L5': param[3] * L_unit, # L5 綁定 L3
        'W7': param[5] * W_unit, 'L7': param[6] * L_unit,
        'W10': param[7] * W_unit, 'L10': param[8] * L_unit
    }
    Cc, Rz, Rf, R_cmfb, I_tail = param[9], param[10], param[11], param[12], param[13]
    R1 = Rf
    
    # ==========================================
    # 1. 節點分類與系統設定
    # ==========================================
    
    # 固定的理想電壓源 (這些節點不參與矩陣求逆)
    fixed_voltages = {
        'VDD': vdd_val, '0': 0.0,
        'inp': vcm_val, 'inn': vcm_val,
        'p1': vcm_val, 'p2': vcm_val, 
        'n1': vcm_val, 'n2': vcm_val
    }

    # 未知電壓節點 (我們要解的變數)
    unknown_nodes = ['1', '2', '3', 's', '4', 'net19', 'net16', 'VBN']
    NUM_VARS = len(unknown_nodes)
    node_to_idx = {n: i for i, n in enumerate(unknown_nodes)}

    # KCL 映射表 (這是通用化的靈魂)
    # 格式: { '物理節點': '對應的方程式列索引(即求解哪個未知數)' }
    kcl_mapping = {
        '1': '1', '2': '2', '3': '3', 's': 's', 
        'net19': 'net19', 'net16': 'net16',
        'VBN': 'VBN',
        # --- 核心魔法 ---
        # p2 和 n2 的電壓已被固定，但我們需要滿足其 KCL。
        # 我們將它們的 KCL 貢獻疊加到 '4' 號方程式，用來求解 V(4)。
        'p2': '4', 
        'n2': '4'  
    }

    # 初始化猜測值 (給予合理的 DC 偏置點)
    v_guess = np.ones(NUM_VARS) * (vdd_val / 2)
    v_guess[node_to_idx['4']] = 0.33  # V_ctrl 初始值
    v_guess[node_to_idx['s']] = 0.17  # Vs 初始值
    v_guess[node_to_idx['1']] = vdd_val - 0.29
    v_guess[node_to_idx['2']] = vdd_val - 0.29
    v_guess[node_to_idx['VBN']] = 0.46

    max_iter = 50
    tol = 1e-10

    for it in range(max_iter):
        J = np.zeros((NUM_VARS, NUM_VARS))
        F = np.zeros(NUM_VARS)

        # ------------------------------------------
        # 通用 Stamp 內部函式
        # ------------------------------------------
        def get_v(node):
            """獲取節點當前電壓 (區分固定與未知)"""
            if node in fixed_voltages: return fixed_voltages[node]
            return v_guess[node_to_idx[node]]

        def stamp(physical_node, var_node, val, is_jacobian=True):
            """
            通用的蓋印操作。如果物理節點有映射到方程式，則填入矩陣。
            """
            # 如果這個節點不需要滿足 KCL (例如 VDD, 0, inp 等純電壓源)，則忽略
            if physical_node not in kcl_mapping: return 
            
            # 找到該物理 KCL 實際要貢獻給哪一個方程式 (例如 p2 -> 方程式 4)
            eq_node = kcl_mapping[physical_node]
            row_idx = node_to_idx[eq_node]
            
            if is_jacobian:
                # 只有當擾動的變數是未知數時，才填入雅可比矩陣
                if var_node in node_to_idx:
                    col_idx = node_to_idx[var_node]
                    J[row_idx, col_idx] += val
            else:
                # 填入殘差向量 F (val預期為流出節點的電流)
                F[row_idx] += val

        def add_resistor(n1, n2, r_val):
            g = 1.0 / r_val
            v1, v2 = get_v(n1), get_v(n2)
            i_r = (v1 - v2) * g
            # 殘差
            stamp(n1, None, i_r, False); stamp(n2, None, -i_r, False)
            # 雅可比
            stamp(n1, n1, g, True); stamp(n1, n2, -g, True)
            stamp(n2, n2, g, True); stamp(n2, n1, -g, True)

        def add_current_source(n_out, n_in, i_val):
            """理想電流源 (由 n_out 流向 n_in)"""
            stamp(n_out, None, i_val, False)
            stamp(n_in, None, -i_val, False)

        def add_mosfet(d, g, s, model, w, l, is_pmos=False):
            vd, vg, vs = get_v(d), get_v(g), get_v(s)
            
            if not is_pmos:
                vgs, vds = vg - vs, vd - vs
                id_w = model.lookup_by_vgs('ID_W', l, vgs, vds)
                gm_w = model.lookup_by_vgs('GM', l, vgs, vds) / W_unit
                gds_w = model.lookup_by_vgs('GDS', l, vgs, vds) / W_unit
                ids, gm, gds = id_w * w, gm_w * w, gds_w * w
                
                # NMOS Stamp
                stamp(d, None, ids, False); stamp(s, None, -ids, False)
                stamp(d, d, gds, True); stamp(d, g, gm, True); stamp(d, s, -(gm+gds), True)
                stamp(s, d, -gds, True); stamp(s, g, -gm, True); stamp(s, s, (gm+gds), True)
            else:
                vsg, vsd = vs - vg, vs - vd
                id_w = model.lookup_by_vgs('ID_W', l, vsg, vsd)
                gm_w = model.lookup_by_vgs('GM', l, vsg, vsd) / W_unit
                gds_w = model.lookup_by_vgs('GDS', l, vsg, vsd) / W_unit
                isd, gm, gds = id_w * w, gm_w * w, gds_w * w
                
                # PMOS Stamp (注意電流方向與電壓導數)
                stamp(d, None, -isd, False); stamp(s, None, isd, False)
                stamp(d, d, gds, True); stamp(d, g, gm, True); stamp(d, s, -(gm+gds), True)
                stamp(s, d, -gds, True); stamp(s, g, -gm, True); stamp(s, s, (gm+gds), True)

        # ==========================================
        # 2. 構建電路 (遍歷 Netlist，呼叫 Stamp)
        # ==========================================
        
        # --- 被動元件 ---
        add_resistor('p2', '4', R_cmfb)  # R38
        add_resistor('n2', '4', R_cmfb)  # R41
        add_resistor('p1', 'n2', Rf)     # R34
        add_resistor('n1', 'p2', Rf)     # R35
        add_resistor('1', '3', R_cmfb)   # R39
        add_resistor('2', '3', R_cmfb)   # R40
        add_resistor('inp', 'p1', R1)    # R37
        add_resistor('inn', 'n1', R1)    # R36
        # (這裡省略電容 C0, C1 的並聯電阻 R42, R43，因為 DC 穩態時電容開路，若無直流路徑可忽略，或給予極小導納)

        # --- 第一級 (NMOS 對 + PMOS 負載) ---
        # 關鍵差異1: M10 替換為理想電流源
        #add_current_source('s', '0', I_tail) 
        add_mosfet('s', 'VBN', '0', nch, dims['W10'], dims['L10'], is_pmos=False)
        add_mosfet('VBN', 'VBN', '0', nch, dims['W10'], dims['L10'], is_pmos=False)
        add_current_source('VDD', 'VBN', I_tail)
        
        add_mosfet('1', 'p1', 's', nch_lvt, dims['W1'], dims['L1'], is_pmos=False) # M1
        add_mosfet('2', 'n1', 's', nch_lvt, dims['W1'], dims['L1'], is_pmos=False) # M2
        add_mosfet('1', '3', 'VDD', pch, dims['W3'], dims['L3'], is_pmos=True)     # M3
        add_mosfet('2', '3', 'VDD', pch, dims['W3'], dims['L3'], is_pmos=True)     # M4

        # --- 第二級 ---
        add_mosfet('p2', '1', 'VDD', pch, dims['W5'], dims['L5'], is_pmos=True)   # M5
        add_mosfet('n2', '2', 'VDD', pch, dims['W5'], dims['L5'], is_pmos=True)   # M6
        add_mosfet('p2', '4', '0', nch, dims['W7'], dims['L7'], is_pmos=False)     # M7
        add_mosfet('n2', '4', '0', nch, dims['W7'], dims['L7'], is_pmos=False)     # M8

        # 注意：我們沒有加入 I35(4, 0)。因為方程式 4 已經被我們挪用去滿足 p2/n2 的 KCL 了。

        # ==========================================
        # 3. 矩陣求解 (Newton-Raphson)
        # ==========================================
        res_norm = np.linalg.norm(F)
        if res_norm < tol:
            break
            
        try:
            delta_v = np.linalg.solve(J + np.eye(NUM_VARS)*1e-15, -F)
            v_guess = v_guess + delta_v * 0.8 # 引入 0.8 的 Damping 避免震盪
        except np.linalg.LinAlgError:
            raise ValueError("矩陣奇異，求解失敗")

    # ==========================================
    # 4. 後處理 (收回丟失的信息)
    # ==========================================
    final_v = {n: v_guess[i] for n, i in node_to_idx.items()}
    final_v.update(fixed_voltages)

    # 處理差異1: 反推 M10 的 VBN (VGS)
    # vds10 = final_v['s']
    # id_w_10 = I_tail / (4e-6 * 8) # I_tail / W10_total
    # vgs10 = nch.lookup_vgs_by_idw(id_w_10, L=500e-9, VDS_guess=vds10)
    # final_v['VBN'] = vgs10

    # 處理差異2: 反推 I35 應該是多少
    # 此時所有電壓都已知，我們重新計算流入節點 4 的電流總和，該總和即為需要被 I35 吸走的電流。
    i_R38 = (final_v['p2'] - final_v['4']) / R_cmfb
    i_R41 = (final_v['n2'] - final_v['4']) / R_cmfb
    i_35_needed = i_R38 + i_R41
    
    VGS1 = final_v['inp'] - final_v['s']
    VDS1 = final_v['1'] - final_v['s']
    VGS3 = VDS3 = VGS5 = final_v['VDD'] - final_v['1']
    VDS5 = vdd_val - vcm_val
    VDS7 = vcm_val
    VGS7 = final_v['4']
    VGS10 = final_v['VBN']
    VDS10 = final_v['s']
    i_tail = nch.lookup_by_vgs('ID_W', dims['L10'], VGS10, VDS10) * dims['W10']
    i_2 = pch.lookup_by_vgs('ID_W', dims['L5'], VGS5, VDS5) * dims['W5']
    i_total = i_tail + i_2 * 2
    final_v['vdsat_tail'] = nch.lookup_by_vgs('VDSAT', dims['L10'], VGS10, VDS10)
    op_config = {
        'mosfets':{
            'M1':  {'VGS': VGS1, 'VDS': VDS1, 'W': dims['W1'], 'L': dims['L1'], 'model': nch_lvt, 'type_inl': 'nch'},
            'M2':  {'VGS': VGS1, 'VDS': VDS1, 'W': dims['W1'], 'L': dims['L1'], 'model': nch_lvt, 'type_inl': 'nch'},
            'M3':  {'VGS': VGS3, 'VDS': VDS3, 'W': dims['W3'], 'L': dims['L3'], 'model': pch, 'type_inl': 'pch'},
            'M4':  {'VGS': VGS3, 'VDS': VDS3, 'W': dims['W3'], 'L': dims['L3'], 'model': pch, 'type_inl': 'pch'},
            'M5':  {'VGS': VGS5, 'VDS': VDS5, 'W': dims['W5'], 'L': dims['L5'], 'model': pch, 'type_inl': 'pch'},
            'M6':  {'VGS': VGS5, 'VDS': VDS5, 'W': dims['W5'], 'L': dims['L5'], 'model': pch, 'type_inl': 'pch'},
            'M7':  {'VGS': VGS7, 'VDS': VDS7, 'W': dims['W7'], 'L': dims['L7'], 'model': nch, 'type_inl': 'nch'},
            'M8':  {'VGS': VGS7, 'VDS': VDS7, 'W': dims['W7'], 'L': dims['L7'], 'model': nch, 'type_inl': 'nch'},
            'M10': {'VGS': VGS10, 'VDS': VDS10, 'W': dims['W10'], 'L': dims['L10'], 'model': nch, 'type_inl': 'nch'},
        },
        # 被動元件一併傳遞
        'passives': {'R_f': Rf, 'R_in': R1, 'Rz': Rz, 'Cc': Cc, 'R_cmfb': R_cmfb}
    }
    final_v['currents'] = {'I_tail': i_tail, 'I_2': i_2, 'I_total': i_total, 'I_adj': i_35_needed}
    #pprint(final_v)
    return op_config, final_v

if __name__ == "__main__":
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
    models = {'nch': nch, 'nch_lvt': nch_lvt, 'pch': pch}

    initial_guess = [
        24,       # W1_units
        10,       # L1_units
        24,       # W3_units
        10,       # L3_units
        48,       # W5_units
        12,       # W7_units
        10,       # L7_units
        8,        # W10_units
        50,       # L10_units
        402e-15,  # Cc
        204.6,    # Rz
        4000,     # Rf/R1
        40e3,     # R_cmfb
        78.5e-6  # I_tail
    ]

    op_config, dc = solve_dc_fully_diff(initial_guess, models)
    pprint(op_config)
    pprint(dc)