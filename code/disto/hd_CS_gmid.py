import numpy as np
import pandas as pd
import re
from scipy.interpolate import RegularGridInterpolator

# ================= 參數設定 =================
ID_FILE = 'id.csv'
GM_FILE = 'gm.csv'
GDS_FILE = 'gds.csv'

# 目標操作點與測試條件 (與你先前的 Transient 設定保持一致)
TARGET_VGS = 0.290       # 290 mV
TARGET_VDS = 0.4375      # 437.5 mV
R_LOAD = 10e3            # 10k 負載
V_AMP_IN = 20e-3          # 20mV AC 擺幅
# ============================================

def extract_vgs_from_string(col_name):
    """從字串中穩健地提取 Vgs 的數值"""
    match = re.search(r'Vgs[^\d\.\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name, re.IGNORECASE)
    if match:
        return float(match.group(1))
    match = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name)
    if match:
        return float(match.group(1))
    raise ValueError(f"無法從欄位名稱解析出 Vgs 數值: {col_name}")

def load_cadence_matrix(filename):
    """讀取 CSV 並提取 2D 矩陣"""
    df = pd.read_csv(filename, sep=None, engine='python')
    
    # vds 為矩陣的列方向 (axis=0)
    vds_array = df.iloc[:, 0].astype(float).values
    
    # vgs 為矩陣的行方向 (axis=1)
    vgs_array = np.array([extract_vgs_from_string(col) for col in df.columns[1:]])
    
    # 數據矩陣 shape: (len(vds), len(vgs))
    data_matrix = df.iloc[:, 1:].astype(float).values
    
    return vds_array, vgs_array, data_matrix

def calculate_stage_distortion():
    print("正在處理 2D LUT 數據...")
    
    vds_arr, vgs_arr, _ = load_cadence_matrix(ID_FILE)
    _, _, gm_matrix = load_cadence_matrix(GM_FILE)
    _, _, gds_matrix = load_cadence_matrix(GDS_FILE)
    
    # 計算步長
    vds_step = vds_arr[1] - vds_arr[0]
    vgs_step = vgs_arr[1] - vgs_arr[0]

    # ================= 核心：計算 Kij 矩陣 =================
    # K10 (gm), K01 (gds) 直接來自模擬器平滑數據
    K10_mat = gm_matrix
    K01_mat = gds_matrix
    
    # 二階項
    K20_mat = 0.5 * np.gradient(K10_mat, vgs_step, axis=1) # gm2
    K02_mat = 0.5 * np.gradient(K01_mat, vds_step, axis=0) # gds2
    K11_mat = np.gradient(K10_mat, vds_step, axis=0)       # gm2d (交叉調變項)
    
    # 三階項
    K30_mat = (1/3) * np.gradient(K20_mat, vgs_step, axis=1) # gm3
    K03_mat = (1/3) * np.gradient(K02_mat, vds_step, axis=0) # gds3
    K21_mat = np.gradient(K20_mat, vds_step, axis=0)         # d(gm2)/dVDS
    K12_mat = np.gradient(K02_mat, vgs_step, axis=1)         # d(gds2)/dVGS
    
    # ================= 提取目標偏壓點的 Kij =================
    points = (vds_arr, vgs_arr)
    target_point = np.array([TARGET_VDS, TARGET_VGS])
    
    try:
        K = {}
        K['K10'] = RegularGridInterpolator(points, K10_mat)(target_point)[0]
        K['K01'] = RegularGridInterpolator(points, K01_mat)(target_point)[0]
        K['K20'] = RegularGridInterpolator(points, K20_mat)(target_point)[0]
        K['K02'] = RegularGridInterpolator(points, K02_mat)(target_point)[0]
        K['K11'] = RegularGridInterpolator(points, K11_mat)(target_point)[0]
        K['K30'] = RegularGridInterpolator(points, K30_mat)(target_point)[0]
        K['K03'] = RegularGridInterpolator(points, K03_mat)(target_point)[0]
        K['K21'] = RegularGridInterpolator(points, K21_mat)(target_point)[0]
        K['K12'] = RegularGridInterpolator(points, K12_mat)(target_point)[0]
    except ValueError as e:
        print(f"\n[錯誤] 插值失敗。請確認目標操作點 ({TARGET_VDS}V, {TARGET_VGS}V) 是否在你匯出的 CSV 範圍內。")
        return

    print(f"\n=== 目標偏壓點 VGS={TARGET_VGS*1000}mV, VDS={TARGET_VDS*1000}mV 本質參數 ===")
    print(f"gm   (K10) : {K['K10']*1000:.4f} mS")
    print(f"gds  (K01) : {K['K01']*1e6:.4f} uS")
    print(f"gm2  (K20) : {K['K20']:.4f} A/V^2")
    print(f"gm3  (K30) : {K['K30']:.4f} A/V^3")
    print(f"gm2d (K11) : {K['K11']:.6f} A/V^2 (交叉項)")

    # ================= 計算 Common-Source 放大級失真 =================
    G_L = 1.0 / R_LOAD
    denominator = G_L + K['K01']
    
    a1 = -K['K10'] / denominator
    numerator_a2 = K['K20'] + K['K11'] * a1 + K['K02'] * (a1**2)
    a2 = -numerator_a2 / denominator
    numerator_a3 = (K['K30'] + K['K21'] * a1 + K['K12'] * (a1**2) + K['K03'] * (a1**3) + 
                    K['K11'] * a2 + 2 * K['K02'] * a1 * a2)
    a3 = -numerator_a3 / denominator

    # 計算諧波
    hd2_lin = 0.5 * np.abs(a2 / a1) * V_AMP_IN
    hd3_lin = 0.25 * np.abs(a3 / a1) * (V_AMP_IN**2)
    
    hd2_dbc = 20 * np.log10(hd2_lin) if hd2_lin > 0 else -999
    hd3_dbc = 20 * np.log10(hd3_lin) if hd3_lin > 0 else -999
    sfdr_dbc = -max(hd2_dbc, hd3_dbc)
    
    print(f"\n=== 單級 CS 放大器靜態失真 (Input={V_AMP_IN*1000}mV, R_L={R_LOAD/1000}k) ===")
    print(f"Gain (a1)  : {a1:.4f} V/V")
    print(f"HD2        : {hd2_dbc:.4f} dBc")
    print(f"HD3        : {hd3_dbc:.4f} dBc")
    print(f"SFDR       : {sfdr_dbc:.4f} dBc")

if __name__ == "__main__":
    calculate_stage_distortion()