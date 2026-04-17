import numpy as np
from scipy.optimize import differential_evolution
from scipy.optimize import root
from lut_engine import MosData
from stb_razavi_para import evaluate_razavi_fully_diff_optimizer_ready 
from hd_razavi_para import evaluate_volterra_sfdr
from noise_razavi_para import evaluate_differential_noise
from pprint import pprint

# =====================================================================
# 1. 系統常數與資料庫初始化
# =====================================================================
VDD = 0.9
V_in_CM = 0.45
V_out_CM = 0.45
#I_TOTAL = 200e-6 # 總功耗預算：200uA

# 載入 LUT 資料庫 (注意 P 管開啟 is_pmos=True)
nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
W_unit = 4e-6
L_unit = 10e-9

# =====================================================================
# 2. 核心約束函數 (將優化變數轉為 MNA 參數)
# =====================================================================
def map_vars_to_params(x, is_debug = False):
    """
    整合了通用 DC 求解器的參數轉換函數
    """
    # [A] 解包幾何尺寸 (注意：x[0], x[2] 等是 Unit 數量)
    dims = {
        'W1': x[0] * W_unit, 'L1': x[1] * L_unit,
        'W3': x[2] * W_unit, 'L3': x[3] * L_unit,
        'W5': x[4] * W_unit, 'L5': x[3] * L_unit, # L5 綁定 L3
        'W7': x[5] * W_unit, 'L7': x[6] * L_unit,
        'W10': x[7] * W_unit, 'L10': x[8] * L_unit
    }
    Cc, Rz, Rf, R_cmfb, I_tail = x[9], x[10], x[11], x[12], x[13]
    R1 = Rf

    # [B] 調用通用 DC 求解器 (這裡定義內部函數以方便存取外部變數)
    def kcl_system(v_guess):
        v1_dc, v_s, v_ctrl = v_guess
        # 護欄
        v1_dc = np.clip(v1_dc, 0.05, VDD - 0.05)
        v_s   = np.clip(v_s, 0.01, V_in_CM - 0.05)
        v_ctrl = np.clip(v_ctrl, 0.05, VDD - 0.05)
        # M3 (Diode connected PMOS 第一級負載)
        vgs3 = VDD - v1_dc
        vds3 = vgs3
        # M1 (NMOS 輸入對)
        vgs1 = V_in_CM - v_s
        vds1 = v1_dc - v_s
        # M5 (第二級 PMOS 輸入)
        vgs5 = VDD - v1_dc  # 與 M3 共用 Gate
        vds5 = VDD - V_out_CM
        # M7 (第二級 NMOS 負載)
        vgs7 = v_ctrl
        vds7 = V_out_CM
        # 計算電流 (與之前的 dc_solver 邏輯一致)
        id3 = pch.lookup_by_vgs('ID_W', dims['L3'], vgs3, vds3) * dims['W3']
        id1 = nch_lvt.lookup_by_vgs('ID_W', dims['L1'], vgs1, vds1) * dims['W1']
        id5 = pch.lookup_by_vgs('ID_W', dims['L5'], vgs5, vds5) * dims['W5']
        id7 = nch.lookup_by_vgs('ID_W', dims['L7'], vgs7, vds7) * dims['W7']
        
        err_node1 = id3 - id1
        err_node_s = (2 * id1) - I_tail
        err_node_out = id5 - id7 - (V_out_CM - v_ctrl)/R_cmfb
        return [err_node1, err_node_s, err_node_out]
    
    # 執行求解
    # 給定合理的初始猜測值 [V_1_DC, V_s, V_ctrl]
    initial_guess = [VDD - 0.29, 0.17,  0.33]
    sol = root(kcl_system, initial_guess) #, method='hybr'

    if not sol.success:
        # 如果 DC 不收斂，直接拋出錯誤，讓 objective 捕獲並給予懲罰
        raise ValueError("DC Convergence Failed")

    v1_dc, v_s, v_ctrl = sol.x

    # [C] 根據收斂後的電壓，提取最終 AC 參數
    # 計算各管子的偏置狀態
    # M1
    VGS1, VDS1 = V_in_CM - v_s, v1_dc - v_s
    # M3
    VGS3, VDS3 = VDD - v1_dc, VDD - v1_dc
    # M5
    VGS5, VDS5 = VDD - v1_dc, VDD - V_out_CM
    # M7
    VGS7, VDS7 = v_ctrl, V_out_CM
    # M10 (假設 Gate 電壓由 I_tail 決定，這步可再精化)
    VDS10 = v_s
    VGS10 = nch.lookup_vgs_by_idw(I_tail/dims['W10'], dims['L10'], VDS10)
    
    # 提取所有 AC 參數 (以 W=W_unit 提取後乘以 Unit 數量)
    # 這裡的邏輯是: LUT 的參數是基於 W=4u 算出來的真實值 (不是密度) 
    # 所以 (lookup_value / W_ref) * W_actual = lookup_value * x[n]
    
    def get_ac(lut, W, L, VGS, VDS):
        return {
            'gm': lut.lookup_by_vgs('GM', L, VGS, VDS) * W / W_unit,
            'gds': lut.lookup_by_vgs('GDS', L, VGS, VDS) * W / W_unit,
            'Cgs': lut.lookup_by_vgs('CGS_W', L, VGS, VDS) * W,
            'Cgd': lut.lookup_by_vgs('CGD_W', L, VGS, VDS) * W,
            'Cdb': lut.lookup_by_vgs('CDB_W', L, VGS, VDS) * W,
            'gmbs': lut.lookup_by_vgs('GMBS', L, VGS, VDS) * W if 'GMBS' in lut.grids else 0
        }
    
    p1 = get_ac(nch_lvt, dims['W1'], dims['L1'], VGS1, VDS1)
    p3 = get_ac(pch, dims['W3'], dims['L3'], VGS3, VDS3)
    p5 = get_ac(pch, dims['W5'], dims['L5'], VGS3, VDS5) # M5 的 VGS 就是 VGS3
    p7 = get_ac(nch, dims['W7'], dims['L7'], VGS7, VDS7) # M7 的 VGS 會由 CMFB 決定，這裡以 V_out_CM 近似
    p10 = get_ac(nch, dims['W10'], dims['L10'], VGS10, VDS10)

    actual_I_total = I_tail + 2 * (pch.lookup_by_vgs('ID_W', dims['L5'], VGS5, VDS5) * dims['W5'])

    test_params = {
        'gm1': p1['gm'], 'gds1': p1['gds'], 'gmbs1': p1['gmbs'],
        'Cgs1': p1['Cgs'], 'Cgd1': p1['Cgd'], 'Cdb1': p1['Cdb'], 'Csb1': p1['Cdb'],
        'gm3': p3['gm'], 'gds3': p3['gds'], 'Cgs3': p3['Cgs'], 'Cgd3': p3['Cgd'], 'Cdb3': p3['Cdb'],
        'gm5': p5['gm'], 'gds5': p5['gds'], 'Cgs5': p5['Cgs'], 'Cgd5': p5['Cgd'], 'Cdb5': p5['Cdb'],
        'gm7': p7['gm'], 'gds7': p7['gds'], 'Cgs7': p7['Cgs'], 'Cgd7': p7['Cgd'], 'Cdb7': p7['Cdb'],
        'gds10': p10['gds'], 'Cdb10': p10['Cdb'],
        'R_f': Rf, 'R_1': R1, 'Rz': Rz, 'Cc': Cc, 'R_cmfb': R_cmfb
    }
    # [C] 封裝為 Volterra 求解器所需的狀態字典
    # 我們將管子物件 (model)、尺寸、偏壓點全部包裝起來
    op_config = {
        'M1':  {'VGS': VGS1, 'VDS': VDS1, 'W': dims['W1'], 'L': dims['L1'], 'model': nch_lvt, 'type_inl': 'nch'},
        'M2':  {'VGS': VGS1, 'VDS': VDS1, 'W': dims['W1'], 'L': dims['L1'], 'model': nch_lvt, 'type_inl': 'nch'},
        'M3':  {'VGS': VGS3, 'VDS': VDS3, 'W': dims['W3'], 'L': dims['L3'], 'model': pch, 'type_inl': 'pch'},
        'M4':  {'VGS': VGS3, 'VDS': VDS3, 'W': dims['W3'], 'L': dims['L3'], 'model': pch, 'type_inl': 'pch'},
        'M5':  {'VGS': VGS5, 'VDS': VDS5, 'W': dims['W5'], 'L': dims['L5'], 'model': pch, 'type_inl': 'pch'},
        'M6':  {'VGS': VGS5, 'VDS': VDS5, 'W': dims['W5'], 'L': dims['L5'], 'model': pch, 'type_inl': 'pch'},
        'M7':  {'VGS': VGS7, 'VDS': VDS7, 'W': dims['W7'], 'L': dims['L7'], 'model': nch, 'type_inl': 'nch'},
        'M8':  {'VGS': VGS7, 'VDS': VDS7, 'W': dims['W7'], 'L': dims['L7'], 'model': nch, 'type_inl': 'nch'},
        'M10': {'VGS': VGS10, 'VDS': VDS10, 'W': dims['W10'], 'L': dims['L10'], 'model': nch, 'type_inl': 'nch'},
        # 被動元件一併傳遞
        'passives': {'R_f': Rf, 'R_in': R1, 'Rz': Rz, 'Cc': Cc, 'R_cmfb': R_cmfb}
    }

    if is_debug:
        #pprint(dims)
        #print(f"Rf/R1: {Rf}")
        #print(f"Rcmfb: {R_cmfb}")
        #print(f"Cc: {Cc}")
        #print(f"Rz: {Rz}")
        id1 = nch_lvt.lookup_by_vgs('ID_W', dims['L1'], VGS1, VDS1) * dims['W1']
        id3 = pch.lookup_by_vgs('ID_W', dims['L3'], VGS3, VDS3) * dims['W3']
        id5 = pch.lookup_by_vgs('ID_W', dims['L5'], VGS5, VDS5) * dims['W5']
        id7 = nch.lookup_by_vgs('ID_W', dims['L7'], VGS7, V_out_CM) * dims['W7']
        print("\n=== DC Operating Point Solved ===")
        print(f"V_1_DC (1st Stage Out) : {v1_dc:.4f} V")
        print(f"V_s    (Tail Node)     : {v_s:.4f} V")
        print(f"VGS7   (CMFB Control)  : {v_ctrl:.4f} V")
        print("-" * 33)
        print(f"ID1    (Input branch)  : {id1*1e6:.2f} uA")
        print(f"ID3    (Input branch)  : {id3*1e6:.2f} uA")
        print(f"ID5    (2nd Stage)     : {id5*1e6:.2f} uA")
        print(f"ID7    (2nd Stage)     : {id7*1e6:.2f} uA")
        print(f"I_tail: {I_tail*1e6} uA")
        #pprint(test_params)

    return test_params, actual_I_total, v_s, op_config

# =====================================================================
# 3. 定義 Cost Function (目標函數)
# =====================================================================
def objective(x):
    try:
        # map_vars_to_params 現在會處理所有 DC 求解
        params, actual_I_total, v_s, op_config = map_vars_to_params(x)
        
        # 如果尾電流管被壓到線性區
        if v_s < 0.12:
            return 1e12

        # 進行 AC 評估
        f_arr = np.logspace(4, 10, 200)
        _, _, _, dc_gain, ugf, pm, gm = evaluate_razavi_fully_diff_optimizer_ready(params, f_arr)
        metrics = evaluate_volterra_sfdr(op_config, V_AMP_IN = 0.316, FIN=53/512*100e6)
        
    except ValueError:
        # 如果 DC 不收斂，回傳極大 Cost，讓優化器避開這組尺寸
        return 1e12
    except Exception:
        return 1e12

    # 計算 Cost (加上你強化的懲罰邏輯)
    # cost = 1000 * (x[9] / 1e12) # Cc
    cost = 0
    if pm < 60:
        cost += 1e12 * (60 - pm)**2
    
    if gm < 15:
        cost += 1e12 * (15 - pm)**2
    
    hd3 = np.abs(metrics['HD3_dBc'])
    if hd3 < 75:
        cost += 500 * (75 - hd3)**2
    
    # 2. DC 增益要求 > 45 dB
    # if dc_gain < 45:
    #     cost += 1000 * ((45 - dc_gain) / 45)**2
        
    #3. UGF 頻寬要求 > 350 MHz
    # if ugf < 350e6:
    #     cost += 2000 * ((350e6 - ugf) / 350e6)**2
    
    if actual_I_total > 200e-6:
        cost += 500 * ((actual_I_total - 200e-6) / 200e-6)**2

    #print(f"PM: {pm}, GM: {gm}, Gain: {dc_gain}, GBW: {ugf}, cost: {cost}")
    return cost

# =====================================================================
# 4. 執行全局優化 (Differential Evolution)
# =====================================================================
if __name__ == "__main__":
    print("Initializing Global Optimizer...")
    
    # 設定搜索邊界 (Bounds)
    # Unit 數量: 1 到 20 (即 W 從 4u 到 80u)
    # L 長度: 從你的 LUT 範圍選擇，例如 100n 到 500n
    # Cc: 40fF 到 5pF
    # Rz: 100 歐姆 到 10k 歐姆
    # Rf, R1: 1k 歐姆 到 10k 歐姆
    
    bounds = [
        (1, 20),           # W1_units
        (3, 50),           # L1_units
        (1, 20),           # W3_units
        (3, 90),           # L3_units
        (1, 40),           # W5_units
        (1, 40),           # W7_units
        (3, 90),           # L7_units
        (2, 20),           # W10_units
        (20, 50),          # L10_units
        (80e-15, 2e-12),   # Cc
        (10, 10000),       # Rz
        (1000, 10000),     # Rf/R1
        (40e3, 100e3),     # R_cmfb
        (20e-6, 100e-6)    # I_tail
    ]

    # 定義整數限制陣列 (對應你的 14 個變數)
    # x[0]:W1, x[2]:W3, x[4]:W5, x[5]:W7, x[7]:W10 是 Unit 數量，設為 True
    # 其他如 L, Cc, Rz, Rf, R1, I_tail 是連續物理值，設為 False
    integrality_mask = [
        True,   # x[0] W1_units
        True,   # x[1] L1_units
        True,   # x[2] W3_units
        True,   # x[3] L3_units
        True,   # x[4] W5_units
        True,   # x[5] W7_units
        True,   # x[6] L7_units
        True,   # x[7] W10_units
        True,   # x[8] L10_units
        False,  # x[9] Cc
        False,  # x[10] Rz
        False,  # x[11] Rf/R1
        False,  # x[12] R_cmfb
        False   # x[13] I_tail
    ]

    # 1. 定義你認為合理的「專家經驗值」 (Initial Guess)
    # 注意順序必須與 bounds 一致 (x[0]~x[13])
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
    # initial_guess = [
    #     4,        # W1_units
    #     6,        # L1_units
    #     9,        # W3_units
    #     68,       # L3_units
    #     27,       # W5_units
    #     15,       # W7_units
    #     35,       # L7_units
    #     6,        # W10_units
    #     24,       # L10_units
    #     170e-15,   # Cc
    #     975,   # Rz
    #     6416.1,   # Rf/R1
    #     69286.8,  # R_cmfb
    #     46.6e-6   # I_tail
    # ]

    para_length = len(bounds)

    # 2. 建立初始群體矩陣
    # DE 的群體大小預設是 popsize * len(bounds)
    pop_size_total = 15 * para_length  # popsize=15, 變數=13
    init_pop = np.zeros((pop_size_total, para_length))

    # 將第一行設為你的初始猜測
    init_pop[0] = initial_guess
    # 剩下的行用隨機填滿 (在 bounds 範圍內)
    for i in range(1, pop_size_total):
        for j in range(para_length):
            low, high = bounds[j]
            init_pop[i, j] = np.random.uniform(low, high)

    # 3. 執行優化時傳入 init
    # result = differential_evolution(objective, bounds, init=init_pop, integrality=integrality_mask, popsize=15, maxiter=50, disp=True)

    # print("\n" + "="*50)
    # print("🎉 Optimization Complete! 🎉")
    # print("="*50)
    
    # # 解包最優解
    # best_x = result.x
    # print(f"Optimal W1   : {best_x[0] * 4} um  | L1  : {best_x[1]*10:.0f} nm")
    # print(f"Optimal W3   : {best_x[2] * 4} um  | L3  : {best_x[3]*10:.0f} nm")
    # print(f"Optimal W5   : {best_x[4] * 4} um  | L5  : {best_x[3]*10:.0f} nm (Bound to L3)")
    # print(f"Optimal W7   : {best_x[5] * 4} um  | L7  : {best_x[6]*10:.0f} nm")
    # print(f"Optimal W10  : {best_x[7] * 4} um  | L10 : {best_x[8]*10:.0f} nm")
    # print(f"Optimal Cc   : {best_x[9]*1e12:.2f} pF")
    # print(f"Optimal Rz   : {best_x[10]:.1f} Ohm")
    # print(f"Optimal Rf/R1: {best_x[11]:.1f} Ohm")
    # print(f"Optimal Rcmfb: {best_x[12]:.1f} Ohm")
    # print(f"Optimal I_tail : {best_x[13]*1e6:.1f} uA")

    # pprint(best_x)
    final_params, final_I_total, tail_vds, op_config = map_vars_to_params(initial_guess, is_debug=True)
    # pprint(op_config)
    metrics = evaluate_volterra_sfdr(op_config, V_AMP_IN = 0.316, FIN=53/512*100e6)
    f_arr = np.logspace(4, 10, 500)
    # pprint(final_params)
    _, g_db, p_deg, final_dc, final_ugf, final_pm, final_gm = evaluate_razavi_fully_diff_optimizer_ready(final_params, f_arr, is_plot=False)
    noise = evaluate_differential_noise(op_config, [10e3])
    
    
    print("\n--- Verified AC Performance ---")
    print(f"DC Gain     : {final_dc:.2f} dB")
    print(f"Phase Margin: {final_pm:.2f} Deg")
    print(f"Gain Margin : {final_gm:.2f} dB")
    print(f"UGF         : {final_ugf/1e6:.2f} MHz")
    print(f"Tail VDS    : {tail_vds:.2f} V")
    print(f"Total Current : {final_I_total*1e6:.2f} uA") # 這裡順便把功耗印出來
    print(f"=== Razavi 全差分 2-Stage OTA 閉環靜態失真 ===")
    print(f"閉環差分增益 (A1): {metrics['Gain_V_V']:.5f} V/V")
    print(f"HD2        : {metrics['HD2_dBc']:.2f} dBc (完美對稱抵消)")
    print(f"HD3        : {metrics['HD3_dBc']:.2f} dBc")
    print(f"SFDR       : {-metrics['HD3_dBc']:.2f} dBc")
    print(f"OIP3       : {metrics['OIP3_dBm']:.2f} dBm")
    print(f"=== Razavi 全差分 2-Stage OTA noise ===")
    print(f"out noise : {noise['S_out_diff'][0]} ")