import numpy as np
from scipy.optimize import differential_evolution, NonlinearConstraint
from functools import lru_cache
from scipy.optimize import root
from lut_engine import MosData
from circuit_netlist import CircuitNetlist
from utils import load_circuit_config
from dc_analysis import solve_dc
from stb_analysis import stb_analysis 
from hd_analysis import hd_analysis
from noise_analysis import noise_analysis
from pprint import pprint

# =====================================================================
# 2. 緩存模擬引擎 (避免重複計算)
# =====================================================================
# 將 maxsize 設為 1，這意味著它只緩存「當前這一個候選解」的結果。
# 當 objective 執行完後，constraints_func 調用相同的 x 時，會直接命中緩存。 (此處應為 constraints_func 先執行，objective 後執行)
# 當優化器換到下一個樣本點時，舊緩存會被拋棄。
# 這樣既消除了重複計算，又將記憶體佔用降到了最低（僅存儲一個解的結果）。
@lru_cache(maxsize=1)
def run_full_simulation(x_tuple):
    """
    執行所有昂貴的模擬計算並緩存結果。
    x_tuple: 將 numpy array 轉為 tuple 以利緩存。
    """
    x = np.array(x_tuple)
    try:
        op_config, dc = solve_dc(netlist, models, config['DC'])
        stb = stb_analysis(op_config, config['STB'], is_plot=True)
        hd = hd_analysis(op_config, config['hd'], V_AMP_IN = 0.316, FIN=53/512*100e6)
        noise = noise_analysis(op_config, config['noise'], [1e6, 10e6])
        
        return {
            'success': True,
            'op_config': op_config,
            'dc': dc,
            'stb': stb,
            'hd': hd,
            'noise': noise
        }
    except Exception:
        return {'success': False}

# =====================================================================
# 2. 定義約束函數 (Hard Constraints)
# =====================================================================
def constraints_func(x):
    """
    定義必須滿足的硬性指標：tail_vds > 0.16, pm > 60, gm > 15
    """
    res = run_full_simulation(tuple(x))
    if not res['success']:
        # 回傳值需與 lb=[0.16, 60, 15] 保持合理的「懲罰距離」
        # -1.0 對應 VDS (量級 0.1V)
        # -180.0 對應 PM (量級 10~100 deg)
        # -100.0 對應 GM (量級 10~100 dB)
        return [-1.0, -180.0, -100.0] 

    dc, stb, op_config = res['dc'], res['stb'], res['op_config']
    # 新的 tail_vds 約束：dc['s'] >= dc['tail_vdsat']
    # 等價於 dc['s'] - dc['tail_vdsat'] >= 0
    tail_vds_violation = dc['s'] - op_config['saved_op']['M10']['VDSAT']

    return [
        tail_vds_violation, # tail_vds - tail_vdsat
        stb['pm'],      # Phase Margin
        stb['gm']       # Gain Margin
    ]

# 定義約束對象：lb (下限), ub (上限)
# tail_vds - tail_vdsat >= 0, pm > 60, gm > 15
nlc = NonlinearConstraint(constraints_func, lb=[0.0, 60, 15], ub=[np.inf, np.inf, np.inf])

# =====================================================================
# 3. 定義 Cost Function (目標函數)
# =====================================================================
def objective(x):
    """
    目標函數只負責「越優越好」的部分。
    採用歸一化加權和，避免某個指標過大主導搜索。
    """
    # 參考基準值 (用於歸一化)
    gain_ref = 40.0
    hd3_ref = 75.0
    noise_ref = 4e-16
    power_ref = 200e-6

    # 權重分配 (可以根據你的優先級調整)
    w_gain  = 1.0  # 增益權重
    w_hd3   = 1.0  # 線性度權重
    w_noise = 1.0  # 雜訊權重
    w_power = 1.0  # 功耗權重

    res = run_full_simulation(tuple(x))
    if not res['success']:
        return 1e12
    
    dc, stb, hd, noise = res['dc'], res['stb'], res['hd'], res['noise']

    i_total = dc['currents']['I_total']
    gain = stb['gain']
    hd3 = np.abs(hd['HD3_dBc'])
    in_noise = noise['irn_psd'][0]

    # 性能 Cost 計算 (越小越好)
    # 加上權重因子 w_...
    cost = - w_gain  * (gain / gain_ref) \
           - w_hd3   * (hd3 / hd3_ref) \
           + w_noise * (in_noise / noise_ref) \
           + w_power * (i_total / power_ref)

    return cost

# =====================================================================
# 4. 執行全局優化 (Differential Evolution)
# =====================================================================
if __name__ == "__main__":
    # print("Initializing Global Optimizer...")
    
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
    # result = differential_evolution(
    #     objective, 
    #     bounds, 
    #     constraints=nlc, 
    #     init=init_pop, 
    #     integrality=integrality_mask,
    #     popsize=30,   
    #     maxiter=150,  # 稍微增加迭代次數，配合更強大的策略
    #     strategy='currenttobest1bin', # 使用平衡型策略，提升全局搜索穩定性
    #     mutation=(0.5, 1),           # 抖動變異率，有助於跳出局部最優
    #     recombination=0.7,           # 保持適當的交叉比例
    #     disp=True,
    #     polish=False  # 關閉局部微調，節省大量時間
    # )

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

    #pprint(best_x)
    
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    nch = MosData('nch_lut.csv', W_ref=4e-6, is_pmos=False)
    pch = MosData('pch_lut.csv', W_ref=4e-6, is_pmos=True)
    models = {'nch': nch, 'nch_lvt': nch_lvt, 'pch': pch}
    # W_unit = 4e-6
    # L_unit = 10e-9
    netlist = CircuitNetlist('buf_razavi.net', dialect='spectre')
    #netlist = CircuitNetlist('buf_razavi_full.mdl', dialect='spectre')
    #netlist = CircuitNetlist('buf_razavi_simp.net', dialect='spectre')

    config = load_circuit_config('config.json')
    op_config, dc = solve_dc(netlist, models, config['DC'])
    stb = stb_analysis(op_config, config['STB'], is_plot=True)
    hd = hd_analysis(op_config, config['hd'], V_AMP_IN = 0.316, FIN=53/512*100e6)
    noise = noise_analysis(op_config, config['noise'], [1e6, 10e6])

    print("\n=== DC ===")
    print(f"Tail VDS      : {dc['s']*1e3:.2f} mV")
    print(f"Tail VDSat    : {op_config['saved_op']['M10']['VDSAT']*1e3:.2f} mV")
    i_m10 = op_config['saved_currents']['M10']; i_m5 = op_config['saved_currents']['M5']
    print(f"Tail Current  : {i_m10*1e6:.2f} uA")
    print(f"Output Current: {i_m5*1e6:.2f} uA")
    print(f"Total Current : {(i_m10 + i_m5*2)*1e6:.2f} uA")
    print(f"CMFB Current  : {op_config['saved_currents']['R38']*1e6:.2f} uA")
    print("=== STB ===")
    print(f"DC Gain     : {stb['dm']['gain']:.2f} dB")
    print(f"Phase Margin: {stb['dm']['pm']:.2f} Deg")
    print(f"Gain Margin : {stb['dm']['gm']:.2f} dB")
    print(f"UGF         : {stb['dm']['ugf']/1e6:.2f} MHz")
    print(f"=== distor ===")
    print(f"閉環差分增益 (A1): {hd['Gain_V_V']:.5f} V/V")
    print(f"HD2        : {hd['HD2_dBc']:.2f} dBc (完美對稱抵消)")
    print(f"HD3        : {hd['HD3_dBc']:.2f} dBc")
    print(f"SFDR       : {hd['SFDR_dBc']:.2f} dBc")
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