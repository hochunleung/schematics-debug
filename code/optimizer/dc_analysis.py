import numpy as np
from lut_engine import MosData
from circuit_netlist import CircuitNetlist
from pprint import pprint
from utils import load_circuit_config

def solve_dc(netlist, models, config_dc):
    vdd_val, vcm_val, is_include_stb = config_dc['vdd_val'], config_dc['vcm_val'], config_dc['is_include_stb']
    dict_models = {
        'nch_mac': models['nch'],
        'pch_mac': models['pch'],
        'nch_lvt_mac': models['nch_lvt']
    }

    # 固定的理想電壓源 (這些節點不參與矩陣求逆)
    fixed_voltages = config_dc['fixed_voltages']
    ignored_isources = config_dc['solver_settings']['ignored_isources']
    # 未知電壓節點 (我們要解的變數)
    all_nodes = netlist.nodes
    unknown_nodes = sorted(list(set(all_nodes) - set(fixed_voltages.keys())))
    if is_include_stb: unknown_nodes.extend(['I_Vprb', 'I_Vi', 'I_evinj'])
    #unknown_nodes = ['1', '2', '3', 's', '4', 'net19', 'net16', 'VBN']
    NUM_VARS = len(unknown_nodes)
    node_to_idx = {n: i for i, n in enumerate(unknown_nodes)}
    #pprint(unknown_nodes); pprint(node_to_idx)

    # KCL 映射表 (這是通用化的靈魂)
    # 格式: { '物理節點': '對應的方程式列索引(即求解哪個未知數)' }
    kcl_tricks = config_dc['solver_settings']['kcl_tricks']
    kcl_skip = set(kcl_tricks) | set(kcl_tricks.values())
    kcl_mapping = kcl_tricks.copy()
    for n in unknown_nodes:
        if n not in kcl_skip:
            kcl_mapping[n] = n


    # 初始化猜測值 (給予合理的 DC 偏置點)
    v_guess = np.ones(NUM_VARS) * (vdd_val / 2)

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
            
        def add_vsource(n_pos, n_neg, v_dc, i_var):
            """
            n_pos: 正節點 (x)
            n_neg: 負節點 (p)
            v_dc: 0
            i_idx: 分配給這個電壓源電流的矩陣索引

            Netlist: vprb (x p) vsource dc=0 type=dc
            物理意義： Vx - Vp = 0, 且有一個未知電流 I_Vprb 從 x 流向 p
            未知數索引： 假設分配給 I_Vprb的索引是i_idx
            """
            vp = get_v(n_pos)
            vn = get_v(n_neg)
            i_val = get_v(i_var) # 獲取當前疊代的猜測電流
            
            # 1. 對 F 向量的貢獻 (殘差)
            # KCL 貢獻：電流從正端流出，流入負端
            stamp(n_pos, None, i_val, is_jacobian=False)
            stamp(n_neg, None, -i_val, is_jacobian=False)
            # KVL 貢獻：方程式殘差為 (V+ - V-) - V_dc
            stamp(i_var, None, (vp - vn) - v_dc, is_jacobian=False)
            
            # 2. 對 J 矩陣的貢獻 (偏導數)
            # 對 KCL 方程式求 d(F)/d(I_idx)
            stamp(n_pos, i_var, 1.0, is_jacobian=True)
            stamp(n_neg, i_var, -1.0, is_jacobian=True)
            # 對 KVL 方程式求 d(F)/d(V)
            stamp(i_var, n_pos, 1.0, is_jacobian=True)
            stamp(i_var, n_neg, -1.0, is_jacobian=True)

        def add_current_source(n_out, n_in, i_val):
            """理想電流源 (由 n_out 流向 n_in)"""
            stamp(n_out, None, i_val, False)
            stamp(n_in, None, -i_val, False)

        def add_vcvs(n_pos, n_neg, ctrl_pos, ctrl_neg, gain, i_var):
            """
            Netlist: evinj (n2 n p2 p) vcvs gain=-1
            物理意義： Vn2 - Vn = -1 * (Vp2 - Vp), 且有一個未知電流 I_evinj 從 n2 流向 n
            方程式改寫： (V_n2 - V_n) + 1 * (V_p2 - V_p) = 0
            """
            vp = get_v(n_pos)
            vn = get_v(n_neg)
            vcp = get_v(ctrl_pos)
            vcn = get_v(ctrl_neg)
            i_val = get_v(i_var)
            
            # 1. F 向量貢獻
            stamp(n_pos, None, i_val, is_jacobian=False)
            stamp(n_neg, None, -i_val, is_jacobian=False)
            # KVL: (V+ - V-) - gain * (Vc+ - Vc-)
            stamp(i_var, None, (vp - vn) - gain * (vcp - vcn), is_jacobian=False)
            
            # 2. J 矩陣貢獻
            # KCL 部分 (與 vsource 相同)
            stamp(n_pos, i_var, 1.0, is_jacobian=True)
            stamp(n_neg, i_var, -1.0, is_jacobian=True)
            # KVL 部分 (包含控制端的偏導數)
            stamp(i_var, n_pos, 1.0, is_jacobian=True)
            stamp(i_var, n_neg, -1.0, is_jacobian=True)
            stamp(i_var, ctrl_pos, -gain, is_jacobian=True)
            stamp(i_var, ctrl_neg, gain, is_jacobian=True)

        def add_pcccs(n_pos, n_neg, gain, ctrl_i_vars, coeffs):
            """
            n_pos: 0 (流出)
            n_neg: n (流入)
            ctrl_i_vars: [idx_vprb, idx_Vi] (控制電流在矩陣中的索引)
            coeffs: [c0, c1, c2...]

            Netlist: fiinj (0 n) pcccs gain=-1 probes=[ vprb Vi ] coeffs=[ 0 1 1 ]
            輸出電流流入節點 n(從 0 流出）, 大小為： I_out = gain * (c_0 + c_1 I_vprb + c_2 I_Vi)
            物理意義： 這個元件不產生新的 KVL 方程式，它只是單純把其他 Vsource 的電流放大後，注入到特定的節點
            代入數字： I_out = -1 * (0 + 1 * I_vprb + 1 * I_Vi) = -(I_vprb + I_Vi)
            """
            # 計算受控電流大小
            i_out = gain * coeffs[0]
            for k, i_var in enumerate(ctrl_i_vars):
                i_out += gain * coeffs[k+1] * get_v(i_var)
                
            # 1. F 向量 (注入電流)
            stamp(n_pos, None, i_out, is_jacobian=False)
            stamp(n_neg, None, -i_out, is_jacobian=False)
            
            # 2. J 矩陣 (對控制電流的偏微分)
            for k, i_var in enumerate(ctrl_i_vars):
                dIout_dIctrl = gain * coeffs[k+1]
                stamp(n_pos, i_var, dIout_dIctrl, is_jacobian=True)
                stamp(n_neg, i_var, -dIout_dIctrl, is_jacobian=True)
        
        def add_mosfet(d, g, s, model, w, l, is_pmos=False):
            vd, vg, vs = get_v(d), get_v(g), get_v(s)
            
            if not is_pmos:
                vgs, vds = vg - vs, vd - vs
                id_w = model.lookup_by_vgs('ID_W', l, vgs, vds)
                gm_w = model.lookup_by_vgs('GM', l, vgs, vds) / model.W_ref
                gds_w = model.lookup_by_vgs('GDS', l, vgs, vds) / model.W_ref
                ids, gm, gds = id_w * w, gm_w * w, gds_w * w
                
                # NMOS Stamp
                stamp(d, None, ids, False); stamp(s, None, -ids, False)
                stamp(d, d, gds, True); stamp(d, g, gm, True); stamp(d, s, -(gm+gds), True)
                stamp(s, d, -gds, True); stamp(s, g, -gm, True); stamp(s, s, (gm+gds), True)
            else:
                vsg, vsd = vs - vg, vs - vd
                id_w = model.lookup_by_vgs('ID_W', l, vsg, vsd)
                gm_w = model.lookup_by_vgs('GM', l, vsg, vsd) / model.W_ref
                gds_w = model.lookup_by_vgs('GDS', l, vsg, vsd) / model.W_ref
                isd, gm, gds = id_w * w, gm_w * w, gds_w * w
                
                # PMOS Stamp (注意電流方向與電壓導數)
                stamp(d, None, -isd, False); stamp(s, None, isd, False)
                stamp(d, d, gds, True); stamp(d, g, gm, True); stamp(d, s, -(gm+gds), True)
                stamp(s, d, -gds, True); stamp(s, g, -gm, True); stamp(s, s, (gm+gds), True)

        # ==========================================
        # 2. 構建電路 (遍歷 Netlist，呼叫 Stamp)
        # ==========================================
        op_config = {'mosfets': {}, 'R': {}, 'C': {}, 'saved_currents': {}}
        resistors = netlist.get_elements_by_type('R')
        for res in resistors:
            node_p, node_n, val = res['nodes'][0], res['nodes'][1], res['params']['r']
            add_resistor(node_p, node_n, val)
            op_config['R'][res['name']] = {
                'val': val,
                'p': node_p,
                'n': node_n
            }
        
        mosfets = {}
        for mostype in ['NMOS', 'PMOS']:
            mosfets[mostype] = netlist.get_elements_by_type(mostype)
            #pprint(mosfets)
            for m in mosfets[mostype]:
                add_mosfet(m['nodes'][0], m['nodes'][1], m['nodes'][2], dict_models[m['model']], m['params']['w'], m['params']['l'], is_pmos=(mostype == 'PMOS'))
        
        isources = netlist.get_elements_by_type('ISOURCE')
        for isource in isources:
            if isource['name'] in ignored_isources: continue
            #pprint(isource)
            add_current_source(isource['nodes'][0], isource['nodes'][1], isource['params']['dc'])

        if is_include_stb:
            # MNA 探針分支擴充 (Ideal Sources & KCL/KVL Stamps)
            vprb = netlist.get_element_by_name('Vprb')
            vi = netlist.get_element_by_name('Vi')
            evinj = netlist.get_element_by_name('evinj')
            fiinj = netlist.get_element_by_name('fiinj')
            add_vsource(vprb['nodes'][0], vprb['nodes'][1], 0.0, 'I_Vprb')
            add_vsource(vi['nodes'][0], vi['nodes'][1], 0.0, 'I_Vi')
            add_vcvs(evinj['nodes'][0], evinj['nodes'][1], evinj['nodes'][2], evinj['nodes'][3], evinj['params']['gain'], 'I_evinj')
            add_pcccs(fiinj['nodes'][0], fiinj['nodes'][1], fiinj['params']['gain'], ['I_Vprb', 'I_Vi'], [0, 1, 1])
            op_config['stb_probes'] = {
                'Vi': {'p': vi['nodes'][0], 'n': vi['nodes'][1]},
                'Vprb': {'p': vprb['nodes'][0], 'n': vprb['nodes'][1]},
                'evinj': {'p': evinj['nodes'][0], 'n': evinj['nodes'][1], 'ctrl_p': evinj['nodes'][2], 'ctrl_n': evinj['nodes'][3]},
                'fiinj': {'p': fiinj['nodes'][0], 'n': fiinj['nodes'][1]}
            }

        # ==========================================
        # 3. 矩陣求解 (Newton-Raphson)
        # ==========================================
        res_norm = np.linalg.norm(F)
        if res_norm < tol:
            break
            
        try:
            delta_v = np.linalg.solve(J + np.eye(NUM_VARS)*1e-15, -F)
            # --- 新增這行安全帶 ---
            # 限制每步最大變化不超過 0.2，避免查表跑飛
            delta_v = np.clip(delta_v, -0.2, 0.2)
            v_guess = v_guess + delta_v * 0.8 # 引入 0.8 的 Damping 避免震盪
        except np.linalg.LinAlgError:
            raise ValueError("矩陣奇異，求解失敗")

    # ==========================================
    # 4. 後處理 (收回丟失的信息)
    # ==========================================
    final_v = {n: v_guess[i] for n, i in node_to_idx.items()}
    final_v.update(fixed_voltages)

    saved_currents = config_dc['saved_currents']

    for mostype in ['NMOS', 'PMOS']:
        #pprint(mosfets)
        for m in mosfets[mostype]:
            node_d, node_g, node_s, node_b = m['nodes'][0], m['nodes'][1], m['nodes'][2], m['nodes'][3]
            v_d, v_g, v_s = final_v[node_d], final_v[node_g], final_v[node_s]
            VGS = abs(v_g - v_s)
            VDS = abs(v_d - v_s)
            op_config['mosfets'][m['name']] = {
                'VGS': VGS,
                'VDS': VDS,
                'd': node_d,
                'g': node_g,
                's': node_s,
                'b': node_b,
                'W': m['params']['w'],
                'L': m['params']['l'],
                'model': dict_models[m['model']],
                'type_inl': mostype
            }
            if m['name'] in saved_currents:
                ids = dict_models[m['model']].lookup_by_vgs('ID_W', m['params']['l'], VGS, VDS) * m['params']['w']
                op_config['saved_currents'][m['name']] = ids

    for res in resistors:
        if res['name'] in saved_currents:
            node_p, node_n, val = res['nodes'][0], res['nodes'][1], res['params']['r']
            v_p, v_n = final_v[node_p], final_v[node_n]
            op_config['saved_currents'][res['name']] = abs(v_p - v_n) / val
    
    capcitors = netlist.get_elements_by_type('C')
    for cap in capcitors:
        node_p, node_n = cap['nodes'][0], cap['nodes'][1]
        op_config['C'][cap['name']] = {
            'val': cap['params']['c'],
            'p': node_p,
            'n': node_n
        }

    op_config['all_nodes'] = all_nodes
        
    # op_config = {
    #     'mosfets':{
    #         'M1':  {'VGS': VGS1, 'VDS': VDS1, 'W': dims['W1'], 'L': dims['L1'], 'model': nch_lvt, 'type_inl': 'nch'},
    #         'M2':  {'VGS': VGS1, 'VDS': VDS1, 'W': dims['W1'], 'L': dims['L1'], 'model': nch_lvt, 'type_inl': 'nch'},
    #         'M3':  {'VGS': VGS3, 'VDS': VDS3, 'W': dims['W3'], 'L': dims['L3'], 'model': pch, 'type_inl': 'pch'},
    #         'M4':  {'VGS': VGS3, 'VDS': VDS3, 'W': dims['W3'], 'L': dims['L3'], 'model': pch, 'type_inl': 'pch'},
    #         'M5':  {'VGS': VGS5, 'VDS': VDS5, 'W': dims['W5'], 'L': dims['L5'], 'model': pch, 'type_inl': 'pch'},
    #         'M6':  {'VGS': VGS5, 'VDS': VDS5, 'W': dims['W5'], 'L': dims['L5'], 'model': pch, 'type_inl': 'pch'},
    #         'M7':  {'VGS': VGS7, 'VDS': VDS7, 'W': dims['W7'], 'L': dims['L7'], 'model': nch, 'type_inl': 'nch'},
    #         'M8':  {'VGS': VGS7, 'VDS': VDS7, 'W': dims['W7'], 'L': dims['L7'], 'model': nch, 'type_inl': 'nch'},
    #         'M10': {'VGS': VGS10, 'VDS': VDS10, 'W': dims['W10'], 'L': dims['L10'], 'model': nch, 'type_inl': 'nch'},
    #     },
    #     # 被動元件一併傳遞
    #     'passives': {'R_f': Rf, 'R_in': R1, 'Rz': Rz, 'Cc': Cc, 'R_cmfb': R_cmfb}
    # }
    return op_config, final_v

if __name__ == "__main__":
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
    models = {'nch': nch, 'nch_lvt': nch_lvt, 'pch': pch}

    netlist = CircuitNetlist('buf_razavi.net', dialect='spectre')
    #netlist = CircuitNetlist('buf_razavi_full.mdl', dialect='spectre')
    #netlist = CircuitNetlist('buf_razavi_simp.net', dialect='spectre')

    config_dc = load_circuit_config('config.json')
    op_config, dc = solve_dc(netlist, models, config_dc)
    pprint(op_config)
    pprint(dc)
    #pprint(op_config['stb_probes'])