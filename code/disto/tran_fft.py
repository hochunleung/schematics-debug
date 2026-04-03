import numpy as np
import pandas as pd
from scipy.fft import fft
from scipy.interpolate import interp1d

# ================= 參數設定 =================
CSV_FILE = 'tran.csv'
START_TIME = 180e-9        # 180 nS
SAMPLE_RATE = 100e6        # 100 MHz (對應 10ns 步長)
N_SAMPLES = 512            # 取樣點數
M_CYCLES = 53              # 相干取樣的週期數
# ============================================

def calculate_coherent_fft():
    print("正在讀取並處理 Transient 數據...")
    # 讀取 CSV (第一行是 header)
    try:
        df = pd.read_csv(CSV_FILE, header=0)
        t_raw = df.iloc[:, 0].astype(float).values
        # 你的第二欄已經扣除了 DC (VT - VDC)，直接使用
        vout_raw = df.iloc[:, 1].astype(float).values
    except Exception as e:
        print(f"讀取 CSV 失敗: {e}")
        return

    # 1. 建立線性插值函數
    # 這能模擬 Cadence Calculator 中的 sample() 函數，精準抓取特定時間點的值
    interpolator = interp1d(t_raw, vout_raw, kind='linear')

    # 2. 生成 Cadence FFT 所使用的精準時間軸
    t_step = 1.0 / SAMPLE_RATE  # 10 nS
    t_sampled = START_TIME + np.arange(N_SAMPLES) * t_step
    
    # 確認取樣時間沒有超出原始數據範圍
    if t_sampled[-1] > t_raw[-1]:
        print(f"錯誤：取樣結束時間 ({t_sampled[-1]*1e6:.3f} us) 超出原始數據範圍 ({t_raw[-1]*1e6:.3f} us)")
        return

    # 3. 提取 512 個取樣點
    vout_sampled = interpolator(t_sampled)

    # 4. 執行 FFT (相干取樣不需要加窗函數)
    yf = fft(vout_sampled)
    
    # 將 FFT 結果轉換為真實的電壓振幅 (Peak Amplitude)
    # 乘以 2/N，且不使用窗函數補償
    mag = (2.0 / N_SAMPLES) * np.abs(yf)
    
    # 5. 直接透過 Bin 索引獲取諧波 (相干取樣的巨大優勢)
    # Bin 索引 = 諧波次數 * M_CYCLES
    fund_idx = M_CYCLES          # 53
    hd2_idx  = 2 * M_CYCLES      # 106
    hd3_idx  = 3 * M_CYCLES      # 159

    fund_mag = mag[fund_idx]
    hd2_mag  = mag[hd2_idx]
    hd3_mag  = mag[hd3_idx]

    # 轉換為 dBV (Cadence 預設的 dB20)
    fund_db = 20 * np.log10(fund_mag) if fund_mag > 0 else -999
    hd2_db  = 20 * np.log10(hd2_mag) if hd2_mag > 0 else -999
    hd3_db  = 20 * np.log10(hd3_mag) if hd3_mag > 0 else -999

    # 計算相對失真 (dBc)
    hd2_dbc = hd2_db - fund_db
    hd3_dbc = hd3_db - fund_db
    sfdr_dbc = -max(hd2_dbc, hd3_dbc)

    # 輸出結果
    print(f"\n=== 相干 FFT 頻譜分析結果 (N={N_SAMPLES}, fs={SAMPLE_RATE/1e6}MHz) ===")
    print(f"主頻 (Bin {fund_idx}) : {fund_db:.4f} dBV")
    print(f"HD2  (Bin {hd2_idx}): {hd2_db:.4f} dBV")
    print(f"HD3  (Bin {hd3_idx}): {hd3_db:.4f} dBV")

    print(f"\n=== 失真指標 (與 Cadence 頻譜比對) ===")
    print(f"HD2:  {hd2_dbc:.4f} dBc")
    print(f"HD3:  {hd3_dbc:.4f} dBc")
    print(f"SFDR: {sfdr_dbc:.4f} dBc")

if __name__ == "__main__":
    calculate_coherent_fft()