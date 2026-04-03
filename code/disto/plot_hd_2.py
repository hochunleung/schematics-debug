import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ================= 參數設定 =================
CSV_FILE = 'vtc_data.csv'  
CROP_POINTS = 10  # 設定要截去頭尾各幾個資料點 (可以根據需要調整，例如 5 或 20)
# ============================================

def visualize_derivatives():
    # 讀取沒有 Header 的 CSV
    data = pd.read_csv(CSV_FILE, header=None)
    
    # 取第 2 欄為 Vin，第 3 欄為 Vout
    vin = data.iloc[:, 1].astype(float).values
    vout = data.iloc[:, 2].astype(float).values
    
    # 計算步長與微分
    vin_step = vin[1] - vin[0]
    
    d1 = np.gradient(vout, vin_step)
    d2 = np.gradient(d1, vin_step)
    d3 = np.gradient(d2, vin_step)
    
    # --- 加入資料裁切 ---
    # 利用切片 [CROP_POINTS : -CROP_POINTS] 去除陣列頭尾的極端值
    vin_plot = vin[CROP_POINTS:-CROP_POINTS]
    vout_plot = vout[CROP_POINTS:-CROP_POINTS]
    d1_plot = d1[CROP_POINTS:-CROP_POINTS]
    d2_plot = d2[CROP_POINTS:-CROP_POINTS]
    d3_plot = d3[CROP_POINTS:-CROP_POINTS]
    
    # 建立 4 個子圖
    fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    
    axs[0].plot(vin_plot * 1000, vout_plot * 1000, label='Vout (mV)', color='blue')
    axs[0].set_ylabel('Vout (mV)')
    axs[0].legend(loc='upper right')
    axs[0].grid(True, linestyle='--', alpha=0.6)
    
    axs[1].plot(vin_plot * 1000, d1_plot, label='d1 (a1, Gain)', color='orange')
    axs[1].set_ylabel('d1 (V/V)')
    axs[1].legend(loc='upper right')
    axs[1].grid(True, linestyle='--', alpha=0.6)
    
    axs[2].plot(vin_plot * 1000, d2_plot, label='d2 (2*a2)', color='green')
    axs[2].set_ylabel('d2 (1/V)')
    axs[2].legend(loc='upper right')
    axs[2].grid(True, linestyle='--', alpha=0.6)
    
    axs[3].plot(vin_plot * 1000, d3_plot, label='d3 (6*a3)', color='red')
    axs[3].set_ylabel('d3 (1/V^2)')
    axs[3].set_xlabel('Vin (mV)')
    axs[3].legend(loc='upper right')
    axs[3].grid(True, linestyle='--', alpha=0.6)
    
    plt.suptitle("DC Transfer Characteristic & Its Derivatives (Cropped)", fontsize=14)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    visualize_derivatives()