import numpy as np
import pandas as pd

# ================= 參數設定 =================
CSV_FILE = 'vtc_data.csv'  # 你的 Cadence 匯出檔案名稱
V_BIAS = 0.290             # 你的輸入 DC 偏壓點 (290mV)
V_AMP = 0.001              # 你的輸入 AC 振幅 (1mV)
# ============================================

def calculate_distortion():
    # 1. 讀取數據 (假設 CSV 沒有 header，第一列是 Vin，第二列是 Vout)
    # 若有 header，請調整 pd.read_csv 的參數
    data = pd.read_csv(CSV_FILE, header=None, names=['in', 'out'])
    vin = data['in'].values
    vout = data['out'].values

    # 確保 Vin 是均勻分佈的，取得電壓步長
    vin_step = vin[1] - vin[0]

    # 2. 數值微積分計算高階導數
    # 一階導 dVout/dVin (a1)
    d1 = np.gradient(vout, vin_step)
    # 二階導 d^2Vout/dVin^2
    d2 = np.gradient(d1, vin_step)
    # 三階導 d^3Vout/dVin^3
    d3 = np.gradient(d2, vin_step)

    # 3. 找到最接近目標偏壓點 V_BIAS 的索引
    idx = (np.abs(vin - V_BIAS)).argmin()
    
    print(f"=== 操作點資訊 ===")
    print(f"目標偏壓點: {V_BIAS*1000} mV, 實際取樣偏壓點: {vin[idx]*1000:.4f} mV")
    print(f"DC 輸出電壓: {vout[idx]*1000:.4f} mV")

    # 4. 提取該偏壓點的泰勒展開係數
    a1 = d1[idx]
    a2 = d2[idx] / 2.0  # 記得除以 2!
    a3 = d3[idx] / 6.0  # 記得除以 3!

    print(f"\n=== 泰勒多項式係數 ===")
    print(f"a1 (Gain) : {a1:.4f} V/V")
    print(f"a2        : {a2:.4f} 1/V")
    print(f"a3        : {a3:.4f} 1/V^2")

    # 5. 計算失真指標 (線性比例)
    hd2_lin = 0.5 * np.abs(a2 / a1) * V_AMP
    hd3_lin = 0.25 * np.abs(a3 / a1) * (V_AMP ** 2)
    thd_lin = np.sqrt(hd2_lin**2 + hd3_lin**2)

    # 轉換為 dBc
    hd2_dbc = 20 * np.log10(hd2_lin)
    hd3_dbc = 20 * np.log10(hd3_lin)
    thd_db = 20 * np.log10(thd_lin)
    sfdr_dbc = -max(hd2_dbc, hd3_dbc) # SFDR 取最大諧波的反相

    print(f"\n=== 失真指標 (Input Amplitude = {V_AMP*1000} mV) ===")
    print(f"HD2  : {hd2_dbc:.4f} dBc")
    print(f"HD3  : {hd3_dbc:.4f} dBc")
    print(f"THD  : {thd_db:.4f} dB")
    print(f"SFDR : {sfdr_dbc:.4f} dBc")

if __name__ == "__main__":
    calculate_distortion()