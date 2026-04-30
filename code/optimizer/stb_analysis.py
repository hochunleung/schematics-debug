import numpy as np
import matplotlib.pyplot as plt
from lut_engine import MosData
from circuit_netlist import CircuitNetlist
from dc_analysis import solve_dc
from utils import load_circuit_config
from pprint import pprint

def stb_analysis(op_config, stb_config, is_plot=False):
    """
    可參數化的 MNA 求解器，專為 gm/ID LUT 優化器設計。
    :param params: 包含所有小訊號參數與寄生電容的字典 (由 LUT 生成)
    :param f_array: 頻率陣列
    :return: 頻率, Loop Gain (dB), Phase (Deg), UGF, PM
    """
    f_array = np.logspace(4, 10, 500)
    w_array = 2 * np.pi * f_array
    is_run_CM = False
    T_tian = {
        'dm': np.zeros(len(w_array), dtype=complex),
        'cm': np.zeros(len(w_array), dtype=complex) if is_run_CM else None
    }
    
    # =================================================================
    # 1. 矩陣索引映射 (Mapping)
    # 確保 MNA[0] = I_vidm, MNA[1] = V_x，以標準化 Tian's Method
    # =================================================================
    GND = -1 # AC 地的約定索引

    fixed_voltages = stb_config.get('fixed_voltages', [])
    #fixed_voltages = ['VDD', 'VCM', '0', 'inp', 'inn']
    all_nodes = op_config['all_nodes']
    stb_nodes = sorted(list(set(all_nodes) - set(fixed_voltages)))
    stb_nodes.extend(['I_Vprb', 'I_Vi', 'I_evinj'])
    NUM_VARS = len(stb_nodes)
    node_to_idx = {n: i for i, n in enumerate(stb_nodes)}
    #pprint(stb_nodes); pprint(node_to_idx)

    #gmbs1 = 138.7e-6;    #csb1 = 25.3e-15

    # 1. 動態提取 K 字典和寄生電容字典
    K_dict = {}
    C_dict = {}
    for inst, cfg in op_config['mosfets'].items():
        #if inst == 'passives': continue
        m = cfg['model']
        L, W, VGS, VDS = cfg['L'], cfg['W'], cfg['VGS'], cfg['VDS']
        
        # 提取 K 參數 (密度 * 實際寬度)
        K_dict[inst] = {k: m.lookup_by_vgs(k+'_W', L, VGS, VDS) * W for k in ['K10', 'K01', 'K20', 'K02', 'K11', 'K30', 'K03', 'K21', 'K12']}
        if is_plot:
            print(f"[{inst}] gm: {K_dict[inst]['K10']*1000:.3f} mS | gds: {K_dict[inst]['K01']*1e6:.2f} uS")
        
        # 提取寄生電容 (密度 * 實際寬度，並確保 Cdb 為正絕對值)
        C_dict[inst] = {
            'Cgs': m.lookup_by_vgs('CGS_W', L, VGS, VDS) * W,
            'Cgd': m.lookup_by_vgs('CGD_W', L, VGS, VDS) * W,
            'Cdb': np.abs(m.lookup_by_vgs('CDB_W', L, VGS, VDS) * W) 
        }

    for i, w in enumerate(w_array):
        s = 1j * w
        A = np.zeros((NUM_VARS, NUM_VARS), dtype=complex)
        
        # =================================================================
        # 2. 輔助函數 (Stamps) 
        # (閉包函數，直接操作當前頻率點的 A 矩陣)
        # =================================================================
        def add_y(name1, name2, y_val):
            """添加雙端導納 (電阻/電容)"""
            n1 = node_to_idx[name1] if name1 in node_to_idx else GND
            n2 = node_to_idx[name2] if name2 in node_to_idx else GND
            if n1 != GND: A[n1, n1] += y_val
            if n2 != GND: A[n2, n2] += y_val
            if n1 != GND and n2 != GND:
                A[n1, n2] -= y_val
                A[n2, n1] -= y_val

        def add_gm(d_name, g_name, src_name, gm):
            d = node_to_idx[d_name] if d_name in node_to_idx else GND
            g = node_to_idx[g_name] if g_name in node_to_idx else GND
            src = node_to_idx[src_name] if src_name in node_to_idx else GND
            """添加 VCCS (跨導)"""
            if d != GND:
                if g != GND:   A[d, g] += gm
                if src != GND: A[d, src] -= gm
            if src != GND:
                if g != GND:   A[src, g] -= gm
                if src != GND: A[src, src] += gm
                
        def add_gmbs(d_name, src_name, gmbs):
            d = node_to_idx[d_name] if d_name in node_to_idx else GND
            src = node_to_idx[src_name] if src_name in node_to_idx else GND
            """添加體效應 (Source 參考 AC GND)
            數學上等效於一個受 Source 負電壓控制的 gm"""
            add_gm(d, GND, src, gmbs)

        # =================================================================
        # 3. 填入被動元件與寄生網路
        # =================================================================
        # 回授與負載電阻
        for res, nodes in op_config.get('R', {}).items():
            add_y(nodes['p'], nodes['n'], 1/nodes['val'])
        for cap, nodes in op_config.get('C', {}).items():
            add_y(nodes['p'], nodes['n'], s * nodes['val'])
        
        # =================================================================
        # 4. 填入主動元件 (電晶體 Model)
        # =================================================================
        for mos, nodes in op_config['mosfets'].items():
            add_gm(nodes['d'], nodes['g'], nodes['s'], K_dict[mos]['K10'])
            add_y(nodes['d'], nodes['s'], K_dict[mos]['K01'])
            # add_gmbs(nodes['d'], nodes['s'], gmbs1)
            add_y(nodes['g'], nodes['s'], s*C_dict[mos]['Cgs'])
            add_y(nodes['d'], nodes['g'], s*C_dict[mos]['Cgd'])
            add_y(nodes['d'], nodes['b'], s*C_dict[mos]['Cdb'])

        # =================================================================
        # 5. MNA 探針分支擴充 (Ideal Sources & KCL/KVL Stamps)
        # =================================================================
        # Vi_dm 探針電壓源 (I_VIDM 流出 N_X, 進入 N_P2)
        vi = op_config['stb_probes']['Vi']
        vi_p, vi_n = node_to_idx[vi['p']], node_to_idx[vi['n']]
        I_vi = node_to_idx['I_Vi']
        A[vi_p, I_vi] = 1; A[vi_n, I_vi] = -1
        A[I_vi, vi_p] = 1; A[I_vi, vi_n] = -1  # KVL: Vx - Vp2 = Vtest
        
        # vprb 探針電壓源 (I_VPRB 流出 N_X, 進入 N_P)
        vprb = op_config['stb_probes']['Vprb']
        vprb_p, vprb_n = node_to_idx[vprb['p']], node_to_idx[vprb['n']]
        I_vprb = node_to_idx['I_Vprb']  
        A[vprb_p, I_vprb] = 1; A[vprb_n, I_vprb] = -1
        A[I_vprb, vprb_p] = 1; A[I_vprb, vprb_n] = -1  # KVL: Vx - Vp =0
        
        # evinj 壓控電壓源 (I_EVINJ 流出 N_N2, 進入 N_N)
        evinj = op_config['stb_probes']['evinj']
        evinj_p, evinj_n = node_to_idx[evinj['p']], node_to_idx[evinj['n']]
        evinj_ctrl_p, evinj_ctrl_n = node_to_idx[evinj['ctrl_p']], node_to_idx[evinj['ctrl_n']]
        I_evinj = node_to_idx['I_evinj']
        A[evinj_p, I_evinj] = 1; A[evinj_n, I_evinj] = -1
        A[I_evinj, evinj_p] = 1; A[I_evinj, evinj_n] = -1; 
        for mode, CMDM in [('dm', -1), ('cm', 1)] if is_run_CM else [('dm', -1)]:
            A[I_evinj, evinj_ctrl_p] = -CMDM; A[I_evinj, evinj_ctrl_n] = CMDM
            
            # fiinj 電流注入邏輯 (將 I_vprb 與 I_vidm 鏡像到 N_N)
            A[evinj_n, I_vprb] = -CMDM; A[evinj_n, I_vi] = -CMDM

            # =========================================================
            # 6. Tian's Method 雙重交流掃描與提取
            # =========================================================
            # AC1: 電壓注入
            rhs_ac1 = np.zeros(NUM_VARS, dtype=complex); rhs_ac1[I_vi] = 1  
            x1 = np.linalg.solve(A, rhs_ac1)
            ac1_I_Vi = x1[I_vi]; ac1_V_x = x1[vi_p]
            
            # AC2: 電流注入
            rhs_ac2 = np.zeros(NUM_VARS, dtype=complex); rhs_ac2[vi_p] = 1  
            x2 = np.linalg.solve(A, rhs_ac2)
            ac2_I_Vi = x2[I_vi]; ac2_V_x = x2[vi_p]
            
            cross_prod = ac1_I_Vi * ac2_V_x - ac1_V_x * ac2_I_Vi
            T_tian[mode][i] = -1 / (1 - 1 / (2 * cross_prod + ac1_V_x + ac2_I_Vi))

    # =========================================================
    # 7. 後處理與指標回傳 (供優化器使用)
    # =========================================================
    results = {}
    for mode in ['dm', 'cm'] if is_run_CM else ['dm']:
        gain_db = 20 * np.log10(np.abs(T_tian[mode]))
        phase_deg = np.unwrap(np.angle(T_tian[mode], deg=True))
        
        dc_gain = gain_db[0]
        ugf_idx = np.argmin(np.abs(gain_db))
        ugf = f_array[ugf_idx]
        pm = 180 + phase_deg[ugf_idx] 
        gain_margin = 0

        if np.min(phase_deg) <= -180:
            pcf_idx = np.argmin(np.abs(phase_deg + 180))
            gain_margin = -gain_db[pcf_idx]
        
        results[mode] = {
            'gain': dc_gain,
            'pm':  pm,
            'ugf': ugf, 
            'gm': gain_margin
        }

        if is_plot:
            print(f"=== {'Differential' if mode == 'dm' else 'Common'} Mode (Tian's Method) ===")
            print(f"DC Loop Gain: {dc_gain:.2f} dB")
            print(f"Loop UGF:     {f_array[ugf_idx]/1e6:.2f} MHz")
            print(f"Phase Margin: {pm:.2f} Degrees")
            print(f"Gain Margin:  {gain_margin:.2f} dB")

            fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
            ax1.semilogx(f_array, gain_db, 'b', linewidth=2)
            ax1.axhline(0, color='r', linestyle='--')
            ax1.set_ylabel('Loop Gain (dB)')
            ax1.grid(True, which="both", ls="-", alpha=0.5)
            ax1.set_title(f"Standardized True SPICE Mimic (PM = {pm:.1f} Deg)")

            ax2.semilogx(f_array, phase_deg, 'g', linewidth=2)
            ax2.axhline(-180, color='r', linestyle='--')
            ax2.set_ylabel('Phase (Deg)')
            ax2.set_xlabel('Frequency (Hz)')
            ax2.grid(True, which="both", ls="-", alpha=0.5)
            plt.tight_layout()
            plt.show()

    return results

# ================= 測試代碼 =================
if __name__ == "__main__":
    # 載入 LUT 資料庫 (注意 P 管開啟 is_pmos=True)
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
    models = {'nch': nch, 'nch_lvt': nch_lvt, 'pch': pch}

    #netlist = CircuitNetlist('buf_razavi_full.mdl', dialect='spectre')
    netlist = CircuitNetlist('buf_razavi.net', dialect='spectre')

    config = load_circuit_config('config.json')
    op_config, dc = solve_dc(netlist, models, config['DC'])

    stb = stb_analysis(op_config, config['STB'], is_plot=True)

    # print("\n--- Verified AC Performance ---")
    # print(f"DC Gain     : {stb['gain']:.2f} dB")
    # print(f"Phase Margin: {stb['pm']:.2f} Deg")
    # print(f"Gain Margin : {stb['gm']:.2f} dB")
    # print(f"UGF         : {stb['ugf']/1e6:.2f} MHz")