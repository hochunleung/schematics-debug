import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

def evaluate_fully_diff_middlebrook():
    # ==========================================
    # Netlist 完整參數 (含 Cgs1)
    # ==========================================
    gm1 = 1.252e-3; ro1 = 10227.53
    co1 = 20.38e-15 + 20.98e-15 + 18.84e-15 + 49.54e-15 # 109.74fF
    
    gm2 = 2.786e-3; ro2 = 6299
    co2 = 37.48e-15 + 10.09e-15 + 6.685e-15 # 54.255fF
    
    Cc = 804e-15; Rz = 204.6
    cgd1 = 14e-15; cgd2 = 30.13e-15
    Cgs1 = 25.42e-15  # 新增第一級輸入寄生電容
    
    Rf = 4e3; R1 = 4e3
    gmin = 1e-12 
    
    f = np.logspace(3, 10, 5000)
    w = 2 * np.pi * f
    s = 1j * w
    T_mb = np.zeros(len(s), dtype=complex)
    
    used_fallback = False
    
    for i, s_val in enumerate(s):
        y1 = (1/ro1) + s_val * co1
        y2 = (1/ro2) + s_val * co2
        Yc = (s_val * Cc) / (1 + s_val * Rz * Cc)
        Y_cross = Yc + s_val * cgd2
        Y_gd1 = s_val * cgd1
        Y_gs1 = s_val * Cgs1
        Y_Rf = 1 / Rf
        Y_R1 = 1 / R1
        
        # ==========================================
        # 12x12 全差分 Middlebrook MNA 矩陣
        # 變數向量: [Vp1, Vn1, Vp2, Vn2, Vp, Vn, Vx, Vinn, Vinp, I_vidm, I_vprb, I_evinj]
        # 索引座標:   0    1    2    3    4   5   6    7     8      9       10      11
        # ==========================================
        A = np.zeros((12, 12), dtype=complex)
        
        # 0. Node p1 (第一級正輸出) KCL:
        # 節點導納 + 跨接 Y_cross 往 p2 + 跨接 Cgd1 往 inp + G10(依賴 inp, inn)
        A[0, 0] = y1 + Y_cross + Y_gd1
        A[0, 2] = -Y_cross
        A[0, 8] = gm1/2 - Y_gd1  # G10(+)控制端(inp) + C13(往inp)
        A[0, 7] = -gm1/2         # G10(-)控制端(inn)
        
        # 1. Node n1 (第一級負輸出) KCL:
        # 節點導納 + 跨接 Y_cross 往 n2 + 跨接 Cgd1 往 inn + G12(依賴 inp, inn)
        A[1, 1] = y1 + Y_cross + Y_gd1
        A[1, 3] = -Y_cross
        A[1, 8] = -gm1/2         # G12(+)控制端(inp)
        A[1, 7] = gm1/2 - Y_gd1  # G12(-)控制端(inn) + C24(往inn)
        
        # 2. Node p2 (第二級正輸出) KCL:
        # 節點導納 + 跨接 Y_cross 往 p1 + G11(依賴 p1) + 探針電流 I_vidm (流入為負)
        A[2, 2] = y2 + Y_cross
        A[2, 0] = gm2 - Y_cross
        A[2, 9] = -1             
        
        # 3. Node n2 (第二級負輸出) KCL:
        # 節點導納 + 跨接 Y_cross 往 n1 + G13(依賴 n1) + evinj電流 (流出為正)
        A[3, 3] = y2 + Y_cross
        A[3, 1] = gm2 - Y_cross
        A[3, 11] = 1             
        
        # 4. Node p (正向測試探針節點) KCL:
        # Rf(往 inn) + 探針電流 I_vprb (流入為負)
        A[4, 4] = Y_Rf + gmin
        A[4, 7] = -Y_Rf
        A[4, 10] = -1            
        
        # 5. Node n (負向測試探針節點) KCL:
        # Rf(往 inp) + I_evinj (流入) + fiinj電流 (抽出, 依賴 I_vprb, I_vidm)
        A[5, 5] = Y_Rf + gmin
        A[5, 8] = -Y_Rf
        A[5, 11] = -1            
        A[5, 10] = 1; A[5, 9] = 1 # KCL流出: -(-1*(I_vprb + I_vidm)) = I_vprb + I_vidm
        
        # 6. Node x (Middlebrook 測試注入點) KCL:
        # I_vidm(流出) + I_vprb(流出) = I17(注入電流)
        A[6, 9] = 1; A[6, 10] = 1
        
        # 7. Node inn (第一級負輸入) KCL:
        # Rf(往 p) + R1(往地) + Cgs1(往地) + Cgd1(往 n1)
        A[7, 7] = Y_Rf + Y_R1 + Y_gs1 + Y_gd1
        A[7, 4] = -Y_Rf
        A[7, 1] = -Y_gd1
        
        # 8. Node inp (第一級正輸入) KCL:
        # Rf(往 n) + R1(往地) + Cgs1(往地) + Cgd1(往 p1)
        A[8, 8] = Y_Rf + Y_R1 + Y_gs1 + Y_gd1
        A[8, 5] = -Y_Rf
        A[8, 0] = -Y_gd1
        
        # 9. Vi_dm 理想電壓源 KVL (Vx - Vp2 = Vz1)
        A[9, 6] = 1; A[9, 2] = -1
        
        # 10. vprb 理想電壓源 KVL (Vx - Vp = Vz2)
        A[10, 6] = 1; A[10, 4] = -1
        
        # 11. evinj 理想電壓源 KVL (Vn2 - Vn = -1 * (Vp2 - Vp))
        A[11, 2] = 1; A[11, 3] = 1; A[11, 4] = -1; A[11, 5] = -1

        # =========================================================
        # Test 1: Middlebrook 電壓注入測試 (Voltage Test)
        # =========================================================
        rhs_v = np.zeros(12, dtype=complex)
        rhs_v[9] = 1  # 驅動 Vi_dm
        
        x1 = np.linalg.solve(A, rhs_v)
        V_forward = x1[6]  # Vx (驅動前饋端)
        V_return  = x1[2]  # Vp2 (返回端)
        
        Tv = -V_return / V_forward  # 電壓迴路增益
        
        # =========================================================
        # Test 2: Middlebrook 電流注入測試 (Current Test)
        # =========================================================
        rhs_i = np.zeros(12, dtype=complex)
        rhs_i[6] = 1  # 驅動 I17 (注入節點 x)
        
        x2 = np.linalg.solve(A, rhs_i)
        I_return  = x2[9]       # I_vidm (從放大器輸出端返回的電流)
        I_forward = x2[10]      # I_vprb (進入反饋網路的前饋電流)
        
        # =========================================================
        # Middlebrook 雙注入法閉環公式
        # =========================================================
        if np.abs(I_forward) < 1e-30:
            T_mb[i] = Tv
            used_fallback = True
        else:
            Ti = I_return / I_forward  # 電流迴路增益
            T_mb[i] = (Tv * Ti - 1) / (Tv + Ti + 2)

    # =========================================================
    # 指標提取與 Bode Plot 繪製 (對齊 Cadence 極性)
    # =========================================================
    gain_db = 20 * np.log10(np.abs(T_mb))
    phase_deg = np.unwrap(np.angle(T_mb, deg=True))
    
    dc_gain = gain_db[0]
    ugf_idx = np.argmin(np.abs(gain_db))
    ugf_hz = f[ugf_idx]
    
    # 穩定性標準：距 -180 度的餘裕
    pm = 180 + phase_deg[ugf_idx] 

    # ==================== 新增 Gain Margin 提取 ====================
    # 尋找相位穿越 -180 度 (Phase Crossover Frequency, PCF) 的位置
    # 為了避免相位根本沒掉到 -180 度 (例如一階系統) 導致抓錯，加一個判斷
    if np.min(phase_deg) <= -180:
        pcf_idx = np.argmin(np.abs(phase_deg + 180))
        gm = -gain_db[pcf_idx]  # 穩定系統在 PCF 的 Gain 會是負的，GM 通常以正值表示安全餘裕
        gm_str = f"{gm:.2f} dB"
    else:
        gm_str = "Infinite (相位未達 -180 度)"
    # ===============================================================

    print("=== 全差分終極閉環 Op-Amp (Middlebrook 雙注入法 + 完整寄生) ===")
    print(f"DC Loop Gain: {dc_gain:.2f} dB")
    print(f"Loop UGF:     {ugf_hz/1e6:.2f} MHz")
    print(f"Phase Margin: {pm:.2f} Degrees")
    print(f"Gain Margin:  {gm_str}")  # <=== 把 GM 印出來
    print(f"狀態: {'啟用單向 Fallback' if used_fallback else '使用標準 Middlebrook 公式'}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    ax1.semilogx(f, gain_db, 'b', linewidth=2)
    ax1.axhline(0, color='r', linestyle='--')
    ax1.set_ylabel('Loop Gain (dB)')
    ax1.grid(True, which="both", ls="-", alpha=0.5)
    ax1.set_title(f"Middlebrook True Differential Bode Plot (PM = {pm:.1f} Deg)")

    ax2.semilogx(f, phase_deg, 'g', linewidth=2)
    ax2.axhline(-180, color='r', linestyle='--') # 繪製 -180 度危險線
    ax2.set_ylabel('Phase (Deg)')
    ax2.set_xlabel('Frequency (Hz)')
    ax2.grid(True, which="both", ls="-", alpha=0.5)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    evaluate_fully_diff_middlebrook()