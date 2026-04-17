import numpy as np
import matplotlib.pyplot as plt

def evaluate_razavi_fully_diff_optimizer_ready(params, f_array, is_plot=False):
    """
    可參數化的 MNA 求解器，專為 gm/ID LUT 優化器設計。
    :param params: 包含所有小訊號參數與寄生電容的字典 (由 LUT 生成)
    :param f_array: 頻率陣列
    :return: 頻率, Loop Gain (dB), Phase (Deg), UGF, PM
    """
    w_array = 2 * np.pi * f_array
    T_tian = np.zeros(len(w_array), dtype=complex)
    
    # =================================================================
    # 1. 矩陣索引映射 (Mapping)
    # 確保 MNA[0] = I_vidm, MNA[1] = V_x，以標準化 Tian's Method
    # =================================================================
    I_VIDM = 0; N_X = 1; N_P2 = 2; N_P = 3; N_N2 = 4; N_N = 5
    N_P1 = 6; N_N1 = 7; N_1 = 8; N_2 = 9; N_3 = 10; N_4 = 11; N_S = 12
    I_VPRB = 13; I_EVINJ = 14
    
    GND = -1 # AC 地的約定索引

    for i, w in enumerate(w_array):
        s = 1j * w
        A = np.zeros((15, 15), dtype=complex)
        
        # =================================================================
        # 2. 輔助函數 (Stamps) 
        # (閉包函數，直接操作當前頻率點的 A 矩陣)
        # =================================================================
        def add_y(n1, n2, y_val):
            """添加雙端導納 (電阻/電容)"""
            if n1 != GND: A[n1, n1] += y_val
            if n2 != GND: A[n2, n2] += y_val
            if n1 != GND and n2 != GND:
                A[n1, n2] -= y_val
                A[n2, n1] -= y_val

        def add_gm(d, g, src, gm):
            """添加 VCCS (跨導)"""
            if d != GND:
                if g != GND:   A[d, g] += gm
                if src != GND: A[d, src] -= gm
            if src != GND:
                if g != GND:   A[src, g] -= gm
                if src != GND: A[src, src] += gm
                
        def add_gmbs(d, src, gmbs):
            """添加體效應 (Source 參考 AC GND)
            數學上等效於一個受 Source 負電壓控制的 gm"""
            add_gm(d, GND, src, gmbs)

        # =================================================================
        # 3. 填入被動元件與寄生網路
        # =================================================================
        # 回授與負載電阻
        add_y(N_P, N_N1, 1/params['R_f'])
        add_y(N_N, N_P1, 1/params['R_f'])
        add_y(N_P1, GND, 1/params['R_1'])
        add_y(N_N1, GND, 1/params['R_1'])
        
        # CMFB 電阻
        add_y(N_P2, N_4, 1/params['R_cmfb'])
        add_y(N_N2, N_4, 1/params['R_cmfb'])
        add_y(N_1, N_3, 1/params['R_cmfb'])
        add_y(N_2, N_3, 1/params['R_cmfb'])
        
        # 米勒補償分支 (Rz + Cc)
        y_miller = (s * params['Cc']) / (1 + s * params['Rz'] * params['Cc'])
        add_y(N_P2, N_1, y_miller)
        add_y(N_N2, N_2, y_miller)

        # =================================================================
        # 4. 填入主動元件 (電晶體 Model)
        # =================================================================
        # M13, M14 (第一級輸入對)
        add_gm(N_1, N_P1, N_S, params['gm1']); add_gmbs(N_1, N_S, params['gmbs1'])
        add_y(N_1, N_S, params['gds1'])
        add_y(N_P1, N_S, s*params['Cgs1']); add_y(N_P1, N_1, s*params['Cgd1'])
        
        add_gm(N_2, N_N1, N_S, params['gm1']); add_gmbs(N_2, N_S, params['gmbs1'])
        add_y(N_2, N_S, params['gds1'])
        add_y(N_N1, N_S, s*params['Cgs1']); add_y(N_N1, N_2, s*params['Cgd1'])
        add_y(N_S, GND, s*params['Csb1'] * 2) # 源極寄生電容

        # M17, M18 (第一級主動負載)
        add_gm(N_1, N_3, GND, params['gm3'])
        add_y(N_1, GND, params['gds3'] + s*(params['Cdb3'] + params['Cdb1']))
        add_y(N_3, GND, s*params['Cgs3']); add_y(N_3, N_1, s*params['Cgd3'])
        
        add_gm(N_2, N_3, GND, params['gm3'])
        add_y(N_2, GND, params['gds3'] + s*(params['Cdb3'] + params['Cdb1']))
        add_y(N_3, GND, s*params['Cgs3']); add_y(N_3, N_2, s*params['Cgd3'])

        # M15, M16 (第二級輸入)
        add_gm(N_P2, N_1, GND, params['gm5'])
        add_y(N_1, GND, s*params['Cgs5']); add_y(N_P2, N_1, s*params['Cgd5'])
        
        add_gm(N_N2, N_2, GND, params['gm5'])
        add_y(N_2, GND, s*params['Cgs5']); add_y(N_N2, N_2, s*params['Cgd5'])

        # M11, M12 (第二級主動負載)
        add_gm(N_P2, N_4, GND, params['gm7'])
        add_y(N_P2, GND, params['gds5'] + params['gds7'] + s*(params['Cdb5'] + params['Cdb7']))
        add_y(N_4, GND, s*params['Cgs7']); add_y(N_4, N_P2, s*params['Cgd7'])
        
        add_gm(N_N2, N_4, GND, params['gm7'])
        add_y(N_N2, GND, params['gds5'] + params['gds7'] + s*(params['Cdb5'] + params['Cdb7']))
        add_y(N_4, GND, s*params['Cgs7']); add_y(N_4, N_N2, s*params['Cgd7'])

        # M10 (尾電流源)
        add_y(N_S, GND, params['gds10'] + s*params['Cdb10'])

        # =================================================================
        # 5. MNA 探針分支擴充 (Ideal Sources & KCL/KVL Stamps)
        # =================================================================
        # Vi_dm 探針電壓源 (I_VIDM 流出 N_X, 進入 N_P2)
        A[N_X, I_VIDM] = 1; A[N_P2, I_VIDM] = -1
        A[I_VIDM, N_X] = 1; A[I_VIDM, N_P2] = -1 # KVL: Vx - Vp2 = Vtest
        
        # vprb 探針電壓源 (I_VPRB 流出 N_X, 進入 N_P)
        A[N_X, I_VPRB] = 1; A[N_P, I_VPRB] = -1
        A[I_VPRB, N_X] = 1; A[I_VPRB, N_P] = -1  # KVL: Vx - Vp = 0
        
        # evinj 壓控電壓源 (I_EVINJ 流出 N_N2, 進入 N_N)
        A[N_N2, I_EVINJ] = 1; A[N_N, I_EVINJ] = -1
        A[I_EVINJ, N_N2] = 1; A[I_EVINJ, N_N] = -1; A[I_EVINJ, N_P2] = 1; A[I_EVINJ, N_P] = -1
        
        # fiinj 電流注入邏輯 (將 I_vprb 與 I_vidm 鏡像到 N_N)
        A[N_N, I_VPRB] += 1; A[N_N, I_VIDM] += 1

        # =========================================================
        # 6. Tian's Method 雙重交流掃描與提取
        # =========================================================
        # AC1: 電壓注入
        rhs_ac1 = np.zeros(15, dtype=complex); rhs_ac1[I_VIDM] = 1  
        x1 = np.linalg.solve(A, rhs_ac1)
        ac1_I_Vi = x1[I_VIDM]; ac1_V_x = x1[N_X]
        
        # AC2: 電流注入
        rhs_ac2 = np.zeros(15, dtype=complex); rhs_ac2[N_X] = 1  
        x2 = np.linalg.solve(A, rhs_ac2)
        ac2_I_Vi = x2[I_VIDM]; ac2_V_x = x2[N_X]
        
        cross_prod = ac1_I_Vi * ac2_V_x - ac1_V_x * ac2_I_Vi
        T_tian[i] = -1 / (1 - 1 / (2 * cross_prod + ac1_V_x + ac2_I_Vi))

    # =========================================================
    # 7. 後處理與指標回傳 (供優化器使用)
    # =========================================================
    gain_db = 20 * np.log10(np.abs(T_tian))
    phase_deg = np.unwrap(np.angle(T_tian, deg=True))
    
    dc_gain = gain_db[0]
    ugf_idx = np.argmin(np.abs(gain_db))
    ugf = f_array[ugf_idx]
    pm = 180 + phase_deg[ugf_idx] 
    gain_margin = 0

    if np.min(phase_deg) <= -180:
        pcf_idx = np.argmin(np.abs(phase_deg + 180))
        gain_margin = -gain_db[pcf_idx]
    
    if is_plot:
        print("=== 全差分 Razavi Op-Amp (標準化 Tian's Method) ===")
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


    return f_array, gain_db, phase_deg, dc_gain, ugf, pm, gain_margin

# ================= 測試代碼 =================
if __name__ == "__main__":
    # 這裡的 params 在自動化流程中，將由 gmid LUT 查表後動態生成
    test_params = {
        'gm1': 1.251e-3, 'gds1': 33.83e-6, 'gmbs1': 138.7e-6,
        'Cgs1': 25.38e-15, 'Cgd1': 14e-15, 'Csb1': 25.3e-15, 'Cdb1': 20.38e-15,
        'gm3': 1.201e-3, 'gds3': 39.16e-6,
        'Cgs3': 24.76e-15, 'Cgd3': 18.83e-15, 'Cdb3': 20.98e-15,
        'gm5': 2.784e-3, 'gds5': 81.73e-6,
        'Cgs5': 49.53e-15, 'Cgd5': 30.15e-15, 'Cdb5': 37.48e-15,
        'gm7': 2.401e-3, 'gds7': 52.5e-6,
        'Cgs7': 22.74e-15, 'Cgd7': 12.68e-15, 'Cdb7': 10.09e-15,
        'gds10': 25.32e-6, 'Cdb10': 11.23e-15,
        'R_f': 4e3, 'R_1': 4e3, 'R_cmfb': 40.51e3,
        'Rz': 204.6, 'Cc': 402e-15
    }
    
    f = np.logspace(3, 10, 1000)
    _, gain, phase, dc, ugf, pm, gm = evaluate_razavi_fully_diff_optimizer_ready(test_params, f)
    print(f"DC Gain: {dc:.2f} dB, UGF: {ugf/1e6:.2f} MHz, PM: {pm:.2f} Deg")