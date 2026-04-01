import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

def evaluate_single_ended_middlebrook_optimized():
    # ==========================================
    # Netlist 參數設定 (單邊優化版)
    # ==========================================
    gm1 = 1.252e-3; ro1 = 10227.53
    co1 = (20.38 + 20.98 + 18.84 + 49.54) * 1e-15 # 109.74fF
    
    gm2 = 2.786e-3; ro2 = 6299
    co2 = (37.48 + 10.09 + 6.685) * 1e-15 # 54.255fF
    
    Cc = 804e-15; Rz = 204.6
    cgd1 = 14e-15; cgd2 = 30.13e-15
    
    Rf = 4e3; R1 = 4e3
    gmin = 1e-12 
    
    f = np.logspace(3, 10, 5000)
    w = 2 * np.pi * f
    s = 1j * w
    T_mb = np.zeros(len(s), dtype=complex)
    
    for i, s_val in enumerate(s):
        # 導納計算
        y1 = (1/ro1) + s_val * (co1 + cgd1) # 節點 1 對地總導納
        y2 = (1/ro2) + s_val * co2          # 節點 A 對地總導納
        
        # 降維核心：將 Rz+Cc 支路與 Cgd2 整合成單一跨接導納 Y_cross
        Yc = (s_val * Cc) / (1 + s_val * Rz * Cc)
        Y_cross = Yc + s_val * cgd2
        
        Y_Rf = 1 / Rf
        Y_R1 = 1 / R1
        
        # ==========================================
        # 5x5 Middlebrook MNA 矩陣 (消除 net1)
        # 變數向量: [V1, VA, VB, Vin, I_vi]
        # 索引座標:  0   1   2    3     4
        # ==========================================
        A = np.zeros((5, 5), dtype=complex)
        
        # 0. Node 1 (第一級輸出) KCL:
        # y1*V1 + Y_cross*(V1 - VA) + G2(依賴 in) = 0
        # G2 = gm1 * (0 - Vin) = -gm1 * Vin
        A[0, 0] = y1 + Y_cross
        A[0, 1] = -Y_cross       # 跨接導納往節點 A
        A[0, 3] = -gm1           # G2 受控電流
        
        # 1. Node A (第二級輸出 / 返回端) KCL:
        # y2*VA + Y_cross*(VA - V1) + G3(依賴 1) - I_vi(探針抽出) = 0
        # G3 = gm2 * (V1 - 0)
        A[1, 1] = y2 + Y_cross
        A[1, 0] = gm2 - Y_cross  # G3受控源 + 跨接導納往節點 1
        A[1, 4] = -1             # I_vi 流入 A 節點 (從 B 流過來)
        
        # 2. Node B (測試注入點 / 驅動端) KCL:
        # Y_Rf*(VB - Vin) + I_vi(探針抽出) = I_inj (電流注入源)
        A[2, 2] = Y_Rf
        A[2, 3] = -Y_Rf          # Rf 跨接往節點 in
        A[2, 4] = 1              # I_vi 流出 B 節點
        
        # 3. Node in (原 net2, 輸入回饋節點) KCL:
        # Y_Rf*(Vin - VB) + Y_R1*Vin = 0
        A[3, 3] = Y_Rf + Y_R1
        A[3, 2] = -Y_Rf          # Rf 跨接往節點 B
        
        # 4. Vi 電壓探針 KVL:
        # VB - VA = V_test
        A[4, 2] = 1; A[4, 1] = -1

        # =========================================================
        # Test 1: Middlebrook 電壓注入 (Voltage Test)
        # =========================================================
        rhs_v = np.zeros(5, dtype=complex)
        rhs_v[4] = 1  # 啟動電壓源: 令 V_test = 1V, I_inj = 0
        
        x1 = np.linalg.solve(A, rhs_v)
        V_forward = x1[2]  # VB
        V_return  = x1[1]  # VA
        Tv = -V_return / V_forward
        
        # =========================================================
        # Test 2: Middlebrook 電流注入 (Current Test)
        # =========================================================
        rhs_i = np.zeros(5, dtype=complex)
        rhs_i[2] = 1  # 啟動電流源: 令 I_inj = 1A, V_test = 0
        
        x2 = np.linalg.solve(A, rhs_i)
        I_return  = x2[4]       # I_vi (經由 A 點返回放大器的電流)
        I_forward = 1 - x2[4]   # I_inj - I_vi (流入 Rf 的前饋電流)
        Ti = I_return / I_forward
        
        # =========================================================
        # Middlebrook 雙注入法閉環公式
        # =========================================================
        T_mb[i] = (Tv * Ti - 1) / (Tv + Ti + 2)

    # =========================================================
    # 指標提取與 Bode Plot 繪製
    # =========================================================
    gain_db = 20 * np.log10(np.abs(T_mb))
    phase_deg = np.unwrap(np.angle(T_mb, deg=True))
    
    dc_gain = gain_db[0]
    ugf_idx = np.argmin(np.abs(gain_db))
    pm = 180 + phase_deg[ugf_idx] 

    # Gain Margin (GM) 提取
    if np.min(phase_deg) <= -180:
        pcf_idx = np.argmin(np.abs(phase_deg + 180))
        gm = -gain_db[pcf_idx]
        gm_str = f"{gm:.2f} dB"
    else:
        gm_str = "Infinite (相位未達 -180 度)"

    print("=== 單邊優化版 Op-Amp (Middlebrook 雙注入法) ===")
    print(f"DC Loop Gain: {dc_gain:.2f} dB")
    print(f"Loop UGF:     {f[ugf_idx]/1e6:.2f} MHz")
    print(f"Phase Margin: {pm:.2f} Degrees")
    print(f"Gain Margin:  {gm_str}")

    # 繪圖
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.semilogx(f, gain_db, 'b', linewidth=2)
    ax1.axhline(0, color='r', linestyle='--')
    ax1.set_ylabel('Loop Gain (dB)')
    ax1.grid(True, which="both", ls="-", alpha=0.5)
    ax1.set_title(f"Optimized Single-Ended Bode Plot (PM = {pm:.1f} Deg)")

    ax2.semilogx(f, phase_deg, 'g', linewidth=2)
    ax2.axhline(-180, color='r', linestyle='--')
    ax2.set_ylabel('Phase (Deg)')
    ax2.set_xlabel('Frequency (Hz)')
    ax2.grid(True, which="both", ls="-", alpha=0.5)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    evaluate_single_ended_middlebrook_optimized()