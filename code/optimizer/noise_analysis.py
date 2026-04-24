import numpy as np
from lut_engine import MosData

def noise_analysis(op_config, freqs, T=300):
    """
    通用化的全差分/單端雜訊評估函數。
    :param op_config: 包含電路拓撲、操作點和元件參數的字典。
                      預期包含 'all_nodes', 'mosfets', 'R', 'C', 'series_RC', 'hd_probes'。
    :param freqs: 頻率陣列或列表。
    :param T: 溫度 (Kelvin)。
    :return: 包含頻率、增益、總輸出雜訊、等效輸入雜訊及各元件貢獻的字典。
    """
    # =================================================================
    # 1. 矩陣索引映射 (Mapping)
    # =================================================================
    GND = -1 
    op_config['hd_probes'] = {
        'input': {'p': 'N_INP', 'n': 'N_INN'},
        'output': {'p': 'N_P2', 'n': 'N_N2'}
    }
    all_nodes = op_config['all_nodes']
    fixed_voltages_to_exclude = ['VDD', 'VCM', '0'] 
    
    # 1.1 建立動態節點清單 (不含 AC 地)
    dynamic_nodes = sorted(list(set(all_nodes) - set(fixed_voltages_to_exclude)))
    node_to_idx = {n: i for i, n in enumerate(dynamic_nodes)}
    NUM_NODES = len(dynamic_nodes)

    # 1.2 識別輸入與輸出探針節點
    input_p_node = op_config['hd_probes']['input']['p']
    input_n_node = op_config['hd_probes']['input']['n']
    output_p_node = op_config['hd_probes']['output']['p']
    output_n_node = op_config['hd_probes']['output']['n']

    # 1.3 劃分子矩陣索引 (用於 Partitioning 求解)
    input_nodes = [input_p_node, input_n_node]
    input_idxs = [node_to_idx[n] for n in input_nodes]
    free_nodes = [n for n in dynamic_nodes if n not in input_nodes]
    free_idxs = [node_to_idx[n] for n in free_nodes]

    # =================================================================
    # 2. 提取元件參數
    # =================================================================
    K_dict, C_dict = {}, {}
    for inst, cfg in op_config['mosfets'].items():
        m = cfg['model']
        L, W, VGS, VDS = cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
        K_dict[inst] = {k: m.lookup_by_vgs(k+'_W', L, VGS, VDS) * W for k in ['K10', 'K01']}
        C_dict[inst] = {
            'Cgs': m.lookup_by_vgs('CGS_W', L, VGS, VDS) * W,
            'Cgd': m.lookup_by_vgs('CGD_W', L, VGS, VDS) * W,
            'Cdb': np.abs(m.lookup_by_vgs('CDB_W', L, VGS, VDS) * W) 
        }

    kT4 = 4 * 1.380649e-23 * T

    # 準備結果容器
    results = {
        'freqs': freqs, 'gain': [], 'out_noise_psd': [], 'irn_psd': [],
        'contributions': {comp: [] for comp in list(op_config['mosfets'].keys()) + 
                                              list(op_config.get('R', {}).keys()) + 
                                              list(op_config.get('series_RC', {}).keys())}
    }

    # =================================================================
    # 3. 頻率迭代與矩陣求解
    # =================================================================
    for f in freqs:
        w = 2 * np.pi * f
        s = 1j * w if w != 0 else 1e-15j
        
        # --- A. 建立 Y 矩陣 (MNA) ---
        Y = np.zeros((NUM_NODES, NUM_NODES), dtype=complex)
        
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
            if d != GND:
                if g != GND:   Y[d, g] += gm
                if src != GND: Y[d, src] -= gm
            if src != GND:
                if g != GND:   Y[src, g] -= gm
                if src != GND: Y[src, src] += gm

        # 填入被動元件
        for r_name, r_cfg in op_config.get('R', {}).items():
            add_y(r_cfg['p'], r_cfg['n'], 1/r_cfg['val'])
        for c_name, c_cfg in op_config.get('C', {}).items():
            add_y(c_cfg['p'], c_cfg['n'], s * c_cfg['val'])
        for rc_name, rc_cfg in op_config.get('series_RC', {}).items():
            y_rc = (s * rc_cfg['C']) / (1 + s * rc_cfg['R'] * rc_cfg['C'])
            add_y(rc_cfg['p'], rc_cfg['n'], y_rc)

        # 填入 MOSFET
        for mos_name, mos_cfg in op_config['mosfets'].items():
            d, g, src, b = mos_cfg['d'], mos_cfg['g'], mos_cfg['s'], mos_cfg['b']
            add_gm(d, g, src, K_dict[mos_name]['K10'])
            add_y(d, src, K_dict[mos_name]['K01'])
            add_y(g, src, s * C_dict[mos_name]['Cgs'])
            add_y(d, g, s * C_dict[mos_name]['Cgd'])
            add_y(d, b, s * C_dict[mos_name]['Cdb'])

        # --- B. 建立 SI (雜訊源矩陣) 字典 ---
        SI_dict = {comp: np.zeros((NUM_NODES, NUM_NODES)) for comp in results['contributions'].keys()}
        
        def add_i2(mat, name1, name2, i2_val):
            n1 = node_to_idx[name1] if name1 in node_to_idx else GND
            n2 = node_to_idx[name2] if name2 in node_to_idx else GND
            if n1 != GND: mat[n1, n1] += i2_val
            if n2 != GND: mat[n2, n2] += i2_val
            if n1 != GND and n2 != GND:
                mat[n1, n2] -= i2_val; mat[n2, n1] -= i2_val

        # 填寫 MOSFET 雜訊
        for inst, cfg in op_config['mosfets'].items():
            i2_n = cfg['model'].get_noise_psd(cfg['L'], cfg['VGS'], cfg['VDS'], cfg['W'], f)
            add_i2(SI_dict[inst], cfg['d'], cfg['s'], i2_n)

        # 填寫電阻雜訊
        for r_name, r_cfg in op_config.get('R', {}).items():
            add_i2(SI_dict[r_name], r_cfg['p'], r_cfg['n'], kT4 / r_cfg['val'])

        # 填寫串聯 RC 雜訊 (僅來自 R)
        # for rc_name, rc_cfg in op_config.get('series_RC', {}).items():
        #     y_rc = (s * rc_cfg['C']) / (1 + s * rc_cfg['R'] * rc_cfg['C'])
        #     i2_rc = kT4 * rc_cfg['R'] * (np.abs(y_rc)**2)
        #     add_i2(SI_dict[rc_name], rc_cfg['p'], rc_cfg['n'], i2_rc)

        # --- C. 提取子矩陣與求解 ---
        Y_ff = Y[np.ix_(free_idxs, free_idxs)]
        Y_fi = Y[np.ix_(free_idxs, input_idxs)]
        
        try:
            Z_ff = np.linalg.inv(Y_ff)
        except np.linalg.LinAlgError:
            continue

        # 1. 求解閉環增益
        V_in = np.array([0.5, -0.5]) # 假設差分輸入
        V_free = -Z_ff @ (Y_fi @ V_in)
        
        # 對應到全節點電壓
        V_nodes = np.zeros(NUM_NODES, dtype=complex)
        V_nodes[free_idxs] = V_free
        V_nodes[input_idxs] = V_in
        
        diff_out = V_nodes[node_to_idx[output_p_node]] - V_nodes[node_to_idx[output_n_node]]
        diff_gain = np.abs(diff_out)
        results['gain'].append(diff_gain)

        # 2. 求解雜訊貢獻
        c_out_full = np.zeros(NUM_NODES)
        c_out_full[node_to_idx[output_p_node]] = 1
        c_out_full[node_to_idx[output_n_node]] = -1
        c_out = c_out_full[free_idxs]
        
        H_out = c_out @ Z_ff 

        total_out_psd = 0
        for comp_name, SI_mat in SI_dict.items():
            SI_ff = SI_mat[np.ix_(free_idxs, free_idxs)]
            comp_out_psd = np.real(H_out @ SI_ff @ H_out.conj().T)
            results['contributions'][comp_name].append(comp_out_psd)
            total_out_psd += comp_out_psd
            
        results['out_noise_psd'].append(total_out_psd)
        results['irn_psd'].append(total_out_psd / (diff_gain**2) if diff_gain > 1e-12 else 0)

    return results

if __name__ == "__main__":
    # 測試代碼範例
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)

    test_op_config = {
        'all_nodes': ['N_INP', 'N_INN', 'N_P1', 'N_N1', 'N_S', 'N_1', 'N_2', 'N_3', 'N_P2', 'N_N2', 'N_4', 'mid1', 'mid2', 'VDD', 'VCM', '0'],
        'mosfets': {
            'M1':  {'d': 'N_1', 'g': 'N_P1', 's': 'N_S', 'b': 'N_S', 'VGS': 0.2565, 'VDS': 0.4191, 'W': 96e-6, 'L': 100e-9, 'model': nch_lvt},
            'M2':  {'d': 'N_2', 'g': 'N_N1', 's': 'N_S', 'b': 'N_S', 'VGS': 0.2565, 'VDS': 0.4191, 'W': 96e-6, 'L': 100e-9, 'model': nch_lvt},
            'M3':  {'d': 'N_1', 'g': 'N_3', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.2873, 'W': 96e-6, 'L': 100e-9, 'model': pch},
            'M4':  {'d': 'N_2', 'g': 'N_3', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.2873, 'W': 96e-6, 'L': 100e-9, 'model': pch},
            'M5':  {'d': 'N_P2', 'g': 'N_1', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.4498, 'W': 192e-6,'L': 100e-9, 'model': pch},
            'M6':  {'d': 'N_N2', 'g': 'N_2', 's': 'VDD', 'b': 'VDD', 'VGS': 0.2873, 'VDS': 0.4498, 'W': 192e-6,'L': 100e-9, 'model': pch},
            'M7':  {'d': 'N_P2', 'g': 'N_4', 's': '0', 'b': '0', 'VGS': 0.3301, 'VDS': 0.4502, 'W': 48e-6, 'L': 100e-9, 'model': nch},
            'M8':  {'d': 'N_N2', 'g': 'N_4', 's': '0', 'b': '0', 'VGS': 0.3301, 'VDS': 0.4502, 'W': 48e-6, 'L': 100e-9, 'model': nch},
            'M10': {'d': 'N_S', 'g': 'VCM', 's': '0', 'b': '0', 'VGS': 0.4585, 'VDS': 0.1936, 'W': 32e-6, 'L': 500e-9, 'model': nch},
        },
        'R': {
            'R_f_p': {'p': 'N_P1', 'n': 'N_N2', 'val': 4000},
            'R_f_n': {'p': 'N_N1', 'n': 'N_P2', 'val': 4000},
            'R_in_p': {'p': 'N_P1', 'n': 'N_INP', 'val': 4000},
            'R_in_n': {'p': 'N_N1', 'n': 'N_INN', 'val': 4000},
            'R_cmfb_p2_4': {'p': 'N_P2', 'n': 'N_4', 'val': 40e3},
            'R_cmfb_n2_4': {'p': 'N_N2', 'n': 'N_4', 'val': 40e3},
            'R_cmfb_1_3': {'p': 'N_1', 'n': 'N_3', 'val': 40e3},
            'R_cmfb_2_3': {'p': 'N_2', 'n': 'N_3', 'val': 40e3},
            'Rz1': {'p': 'mid1', 'n': 'N_1', 'val': 204.6},
            'Rz2': {'p': 'mid2', 'n': 'N_2', 'val': 204.6},
        },
        'C': {
            'Cc1': {'p': 'N_P2', 'n': 'mid1', 'val': 402e-15},
            'Cc2': {'p': 'N_N2', 'n': 'mid2', 'val': 402e-15},
        },
        # 'series_RC': {
        #     'miller_p': {'p': 'N_P2', 'n': 'N_1', 'R': 204.6, 'C': 402e-15},
        #     'miller_n': {'p': 'N_N2', 'n': 'N_2', 'R': 204.6, 'C': 402e-15},
        # },
        # 'hd_probes': {
        #     'input': {'p': 'N_INP', 'n': 'N_INN'},
        #     'output': {'p': 'N_P2', 'n': 'N_N2'}
        # }
    }
    
    noise_results = noise_analysis(test_op_config, [10e6])
    print(f"Gain      : {noise_results['gain'][0]:.4f}")
    print(f"Out Noise : {noise_results['out_noise_psd'][0]:.2e} V^2/Hz")
    print(f"In  Noise : {noise_results['irn_psd'][0]:.2e} V^2/Hz")
    
    # 列印前三大雜訊貢獻者
    contributions = {k: v[0] for k, v in noise_results['contributions'].items()}
    sorted_contri = sorted(contributions.items(), key=lambda x: x[1], reverse=True)
    print("\nTop 3 Noise Contributors:")
    for name, val in sorted_contri[:3]:
        print(f"  {name}: {val:.2e} V^2/Hz ({val/noise_results['out_noise_psd'][0]*100:.1f}%)")