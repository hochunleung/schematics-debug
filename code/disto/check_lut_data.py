import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re

# ================= 參數設定 =================
ID_FILE = 'id.csv'
GM_FILE = 'gm.csv'
GDS_FILE = 'gds.csv'
# ============================================

def extract_vgs_from_string(col_name):
    """從字串中穩健地提取 Vgs 的數值，忽略 (A), (S) 等單位"""
    # 尋找 Vgs 後面的浮點數
    match = re.search(r'Vgs[^\d\.\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name, re.IGNORECASE)
    if match:
        return float(match.group(1))
    
    # 如果找不到 Vgs 關鍵字，退而求其次找字串中的第一個數字
    match = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name)
    if match:
        return float(match.group(1))
        
    raise ValueError(f"無法從欄位名稱解析出 Vgs 數值: {col_name}")

def load_cadence_matrix(filename):
    """讀取 CSV 並提取 2D 矩陣"""
    df = pd.read_csv(filename, sep=None, engine='python')
    
    # 1. 提取 VDS 陣列 (第 0 欄)
    vds_array = df.iloc[:, 0].astype(float).values
    
    # 2. 提取 VGS 陣列 (使用正則表達式)
    vgs_array = np.array([extract_vgs_from_string(col) for col in df.columns[1:]])
    
    # 3. 提取數據矩陣
    data_matrix = df.iloc[:, 1:].astype(float).values
    
    return vds_array, vgs_array, data_matrix

def visualize_lut_data():
    print("正在讀取並解析 CSV 數據...")
    
    vds_id, vgs_id, id_mat = load_cadence_matrix(ID_FILE)
    vds_gm, vgs_gm, gm_mat = load_cadence_matrix(GM_FILE)
    vds_gds, vgs_gds, gds_mat = load_cadence_matrix(GDS_FILE)
    
    print(f"成功讀取數據！")
    print(f"VDS 軸長度: {len(vds_id)}, 範圍: {vds_id[0]:.3f}V ~ {vds_id[-1]:.3f}V")
    print(f"VGS 軸長度: {len(vgs_id)}, 範圍: {vgs_id[0]:.3f}V ~ {vgs_id[-1]:.3f}V")
    print("-" * 30)
    
    # 建立 3 個子圖 (Id, gm, gds)
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    
    # 為了避免線條過多造成畫面變成一團黑，我們均勻挑選約 10 條 VGS 曲線來畫
    step = max(1, len(vgs_id) // 10)
    
    for i in range(0, len(vgs_id), step):
        vgs_val = vgs_id[i]
        
        # 繪製 Id (轉換為 uA)
        axs[0].plot(vds_id * 1000, id_mat[:, i] * 1e6, label=f'{vgs_val*1000:.0f}mV')
        # 繪製 gm (轉換為 mS)
        axs[1].plot(vds_gm * 1000, gm_mat[:, i] * 1000, label=f'{vgs_val*1000:.0f}mV')
        # 繪製 gds (轉換為 uS)
        axs[2].plot(vds_gds * 1000, gds_mat[:, i] * 1e6, label=f'{vgs_val*1000:.0f}mV')
        
    # 設定圖表格式
    axs[0].set_title('Id vs VDS')
    axs[0].set_xlabel('VDS (mV)')
    axs[0].set_ylabel('Id (uA)')
    axs[0].grid(True, linestyle='--', alpha=0.6)
    
    axs[1].set_title('gm vs VDS')
    axs[1].set_xlabel('VDS (mV)')
    axs[1].set_ylabel('gm (mS)')
    axs[1].grid(True, linestyle='--', alpha=0.6)
    
    axs[2].set_title('gds vs VDS')
    axs[2].set_xlabel('VDS (mV)')
    axs[2].set_ylabel('gds (uS)')
    axs[2].grid(True, linestyle='--', alpha=0.6)
    
    # 加上共用的圖例
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='center right', title='VGS')
    
    plt.suptitle("Cadence 2D Sweep Data Visualization", fontsize=14)
    plt.tight_layout(rect=[0, 0, 0.92, 1]) # 留出右側空間給圖例
    plt.show()

if __name__ == "__main__":
    visualize_lut_data()