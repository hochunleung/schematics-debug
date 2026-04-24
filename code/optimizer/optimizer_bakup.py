import numpy as np
from scipy.optimize import differential_evolution
from scipy.optimize import root
from lut_engine import MosData
from dc_razavi import solve_dc_fully_diff
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
models = {'nch': nch, 'nch_lvt': nch_lvt, 'pch': pch}
W_unit = 4e-6
L_unit = 10e-9

# =====================================================================
# 3. 定義 Cost Function (目標函數)
# =====================================================================
def objective(x):
    try:
        op_config, dc = solve_dc_fully_diff(x, models)
        stb = evaluate_razavi_fully_diff_optimizer_ready(op_config, is_plot=False)
        hd = evaluate_volterra_sfdr(op_config, V_AMP_IN = 0.316, FIN=53/512*100e6)
        noise = evaluate_differential_noise(op_config, [10e6])
        
    except ValueError:
        # 如果 DC 不收斂，回傳極大 Cost，讓優化器避開這組尺寸
        return 1e12
    except Exception:
        return 1e12
    
    tail_vds = dc['s']
    i_total = dc['currents']['I_total']
    pm, gm, gain, ugf= stb['pm'], stb['gm'], stb['gain'], stb['ugf']
    hd3 = np.abs(hd['HD3_dBc'])
    out_noise = noise['out_noise_psd'][0]
    in_noise = noise['irn_psd'][0]
    # 計算 Cost (加上你強化的懲罰邏輯)
    # cost = 1000 * (x[9] / 1e12) # Cc
    tail_vds_target = 0.16
    pm_target = 60
    gm_target = 15
    hd3_target = 75
    i_total_target = 200e-6
    in_noise_target = 4e-16
    gain_target = 40

    cost = 0
    if tail_vds < tail_vds_target:
        cost += 1e12 * (1 + (tail_vds_target - tail_vds) / tail_vds_target)**2

    if pm < pm_target:
        cost += 1e12 * (1 + (pm_target - pm) / pm_target)**2
    
    if gm < gm_target:
        cost += 1e12 * (1 + (gm_target - gm) / gm_target)**2

    if hd3 < hd3_target:
        cost += 500 * (1 + (hd3_target - hd3) / hd3_target)**2
    
    if gain < gain_target:
        cost += 500 * (1 + (gain_target - gain) / gain_target)**2
    
    if in_noise > in_noise_target:
        cost += 500 * ((in_noise - in_noise_target) / in_noise_target)**2
    
    #3. UGF 頻寬要求 > 350 MHz
    # if ugf < 350e6:
    #     cost += 2000 * ((350e6 - ugf) / 350e6)**2
    
    if i_total > i_total_target:
        cost += 500 * ((i_total - i_total_target) / i_total_target)**2

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
    result = differential_evolution(objective, bounds, init=init_pop, integrality=integrality_mask, popsize=15, maxiter=50, disp=True)

    print("\n" + "="*50)
    print("🎉 Optimization Complete! 🎉")
    print("="*50)
    
    # 解包最優解
    best_x = result.x
    print(f"Optimal W1   : {best_x[0] * 4} um  | L1  : {best_x[1]*10:.0f} nm")
    print(f"Optimal W3   : {best_x[2] * 4} um  | L3  : {best_x[3]*10:.0f} nm")
    print(f"Optimal W5   : {best_x[4] * 4} um  | L5  : {best_x[3]*10:.0f} nm (Bound to L3)")
    print(f"Optimal W7   : {best_x[5] * 4} um  | L7  : {best_x[6]*10:.0f} nm")
    print(f"Optimal W10  : {best_x[7] * 4} um  | L10 : {best_x[8]*10:.0f} nm")
    print(f"Optimal Cc   : {best_x[9]*1e12:.2f} pF")
    print(f"Optimal Rz   : {best_x[10]:.1f} Ohm")
    print(f"Optimal Rf/R1: {best_x[11]:.1f} Ohm")
    print(f"Optimal Rcmfb: {best_x[12]:.1f} Ohm")
    print(f"Optimal I_tail : {best_x[13]*1e6:.1f} uA")

    pprint(best_x)
    op_config, dc = solve_dc_fully_diff(best_x, models)
    stb = evaluate_razavi_fully_diff_optimizer_ready(op_config, is_plot=True)
    hd = evaluate_volterra_sfdr(op_config, V_AMP_IN = 0.316, FIN=53/512*100e6)
    noise = evaluate_differential_noise(op_config, [1e6, 10e6])
    
    print("\n=== DC ===")
    print(f"Tail VDS      : {dc['s']*1e3:.2f} mV")
    print(f"Tail Current  : {dc['currents']['I_tail']*1e6:.2f} uA")
    print(f"Output Current: {dc['currents']['I_2']*1e6:.2f} uA")
    print(f"Total Current : {dc['currents']['I_total']*1e6:.2f} uA")
    print(f"CMFB Current  : {dc['currents']['I_adj']*1e6:.2f} uA")
    print("=== STB ===")
    print(f"DC Gain     : {stb['gain']:.2f} dB")
    print(f"Phase Margin: {stb['pm']:.2f} Deg")
    print(f"Gain Margin : {stb['gm']:.2f} dB")
    print(f"UGF         : {stb['ugf']/1e6:.2f} MHz")
    print(f"=== distor ===")
    print(f"閉環差分增益 (A1): {hd['Gain_V_V']:.5f} V/V")
    print(f"HD2        : {hd['HD2_dBc']:.2f} dBc (完美對稱抵消)")
    print(f"HD3        : {hd['HD3_dBc']:.2f} dBc")
    print(f"SFDR       : {-hd['HD3_dBc']:.2f} dBc")
    print(f"OIP3       : {hd['OIP3_dBm']:.2f} dBm")
    print(f"=== noise ===")
    print(f"@1MHz ------")
    print(f"Gain      : {noise['gain'][0]} ")
    print(f"out noise : {noise['out_noise_psd'][0]} ")
    print(f"in  noise : {noise['irn_psd'][0]} ")
    print(f"@10MHz ------")
    print(f"Gain      : {noise['gain'][1]} ")
    print(f"out noise : {noise['out_noise_psd'][1]} ")
    print(f"in  noise : {noise['irn_psd'][1]} ")