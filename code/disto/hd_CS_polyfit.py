import numpy as np
import pandas as pd

# ================= 參數設定 =================
CSV_FILE = 'vtc_data.csv'  # 確保你的 CSV 已經去除了第一行文字標題
V_BIAS = 0.290             # DC 偏壓點 (290mV)
V_AMP = 0.001              # AC 振幅 (1mV)
FIT_WINDOW = 0.002         # 擬合窗口：取 V_BIAS 上下 2mV 的數據進行擬合
# ============================================

def calculate_distortion_polyfit():
    # 1. 讀取數據 (無 header)
    data = pd.read_csv(CSV_FILE, header=None)

    # 提取 Vin 和 Vout (第 2 欄和第 3 欄)
    vin = data.iloc[:, 1].astype(float).values
    vout = data.iloc[:, 2].astype(float).values

    # 2. 框選擬合區間
    # 泰勒展開是在「工作點附近」的局部逼近，所以我們只取偏壓點附近的數據
    mask = (vin >= V_BIAS - FIT_WINDOW) & (vin <= V_BIAS + FIT_WINDOW)
    vin_window = vin[mask]
    vout_window = vout[mask]

    if len(vin_window) < 4:
        print("錯誤：偏壓點附近的數據點太少，無法進行 3 階擬合。")
        return

    # 3. 座標平移 (轉為小信號模型)
    # 將輸入電壓減去偏壓，讓 x 軸的 0 對應到 V_BIAS
    vin_ss = vin_window - V_BIAS

    # 4. 執行多項式擬合
    # 擬合公式: Vout = a3*(vin_ss)^3 + a2*(vin_ss)^2 + a1*(vin_ss) + a0
    coeffs = np.polyfit(vin_ss, vout_window, 5)

    # np.polyfit 回傳的陣列順序是從最高次項到最低次項
    a3 = coeffs[2]
    a2 = coeffs[3]
    a1 = coeffs[4]
    a0 = coeffs[5]

    print(f"=== 操作點資訊 ===")
    print(f"目標偏壓點: {V_BIAS*1000} mV")
    print(f"DC 輸出電壓 (擬合結果 a0): {a0*1000:.4f} mV")

    print(f"\n=== 泰勒多項式係數 (來自 Polyfit) ===")
    print(f"a1 (Gain) : {a1:.4f} V/V")
    print(f"a2        : {a2:.4f} 1/V")
    print(f"a3        : {a3:.4f} 1/V^2")

    # 5. 計算失真指標
    # 注意：這裡的 a1, a2, a3 已經是泰勒級數展開後的係數，不需要再除以階乘！
    # 公式：HD2 = 0.5 * |a2/a1| * V_AMP, HD3 = 0.25 * |a3/a1| * V_AMP^2
    hd2_lin = 0.5 * np.abs(a2 / a1) * V_AMP
    hd3_lin = 0.25 * np.abs(a3 / a1) * (V_AMP ** 2)
    thd_lin = np.sqrt(hd2_lin**2 + hd3_lin**2)

    # 轉換為 dBc
    hd2_dbc = 20 * np.log10(hd2_lin) if hd2_lin > 0 else -999
    hd3_dbc = 20 * np.log10(hd3_lin) if hd3_lin > 0 else -999
    thd_db = 20 * np.log10(thd_lin) if thd_lin > 0 else -999
    sfdr_dbc = -max(hd2_dbc, hd3_dbc)

    print(f"\n=== 失真指標 (Input Amplitude = {V_AMP*1000} mV) ===")
    print(f"HD2  : {hd2_dbc:.4f} dBc")
    print(f"HD3  : {hd3_dbc:.4f} dBc")
    print(f"THD  : {thd_db:.4f} dB")
    print(f"SFDR : {sfdr_dbc:.4f} dBc")

if __name__ == "__main__":
    calculate_distortion_polyfit()