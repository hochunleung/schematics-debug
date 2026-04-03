import numpy as np
import pandas as pd
import re
from scipy.interpolate import RegularGridInterpolator

# ================= 參數設定 =================
ID_FILE = 'id.csv'
GM_FILE = 'gm.csv'
GDS_FILE = 'gds.csv'
CDD_FILE = 'cdd.csv'
CGD_FILE = 'cgd.csv'
CGG_FILE = 'cgg.csv'

# 目標操作點與測試條件
TARGET_VGS = 0.290       # 290 mV
TARGET_VDS = 0.4375      # 437.5 mV
R_LOAD = 10e3            # 10k 負載
V_AMP_IN = 20e-3          # 20mV AC 擺幅
FIN = 53 / 512 * 100e6   # 輸入頻率 10.3515625 MHz
# ============================================

def extract_vgs_from_string(col_name):
    match = re.search(r'Vgs[^\d\.\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name, re.IGNORECASE)
    if match: return float(match.group(1))
    match = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name)
    if match: return float(match.group(1))
    raise ValueError(f"無法從欄位名稱解析出 Vgs: {col_name}")

def load_cadence_matrix(filename):
    df = pd.read_csv(filename, sep=None, engine='python')
    vds_array = df.iloc[:, 0].astype(float).values
    vgs_array = np.array([extract_vgs_from_string(col) for col in df.columns[1:]])
    data_matrix = df.iloc[:, 1:].astype(float).values
    return vds_array, vgs_array, data_matrix

def calculate_volterra_distortion():
    print(f"正在處理 2D LUT 數據 (包含電容高頻效應, FIN={FIN/1e6:.4f} MHz)...")
    
    # 讀取靜態矩陣
    vds_arr, vgs_arr, _ = load_cadence_matrix(ID_FILE)
    _, _, gm_matrix = load_cadence_matrix(GM_FILE)
    _, _, gds_matrix = load_cadence_matrix(GDS_FILE)
    
    # 讀取動態電容矩陣
    try:
        _, _, cdd_matrix = load_cadence_matrix(CDD_FILE)
        _, _, cgd_matrix = load_cadence_matrix(CGD_FILE)
        # Cadence OP 的電容有時會帶負號 (-dQ/dV)，這裡統一取絕對值表示實體電容大小
        cdd_matrix = np.abs(cdd_matrix)
        cgd_matrix = np.abs(cgd_matrix)
    except Exception as e:
        print(f"讀取電容 CSV 失敗: {e}")
        return
    
    vds_step = vds_arr[1] - vds_arr[0]
    vgs_step = vgs_arr[1] - vgs_arr[0]

    # ================= 1. 計算靜態 Kij =================
    K10_mat = gm_matrix
    K01_mat = gds_matrix
    K20_mat = 0.5 * np.gradient(K10_mat, vgs_step, axis=1)
    K02_mat = 0.5 * np.gradient(K01_mat, vds_step, axis=0)
    K11_mat = np.gradient(K10_mat, vds_step, axis=0)
    K30_mat = (1/3) * np.gradient(K20_mat, vgs_step, axis=1)
    K03_mat = (1/3) * np.gradient(K02_mat, vds_step, axis=0)
    K21_mat = np.gradient(K20_mat, vds_step, axis=0)
    K12_mat = np.gradient(K02_mat, vgs_step, axis=1)

    # ================= 2. 計算動態電容導數 Cij =================
    # 對 VDS 求導，獲得 Cdd 的非線性係數
    Cdd0_mat = cdd_matrix
    Cdd1_mat = np.gradient(Cdd0_mat, vds_step, axis=0)
    Cdd2_mat = 0.5 * np.gradient(Cdd1_mat, vds_step, axis=0)
    Cgd0_mat = cgd_matrix # Cgd 的非線性在這裡相對 Cdd 較小，暫且只取線性值

    # ================= 3. 提取目標偏壓點 =================
    points = (vds_arr, vgs_arr)
    target_point = np.array([TARGET_VDS, TARGET_VGS])
    
    def extract(mat):
        return RegularGridInterpolator(points, mat)(target_point)[0]

    K = {
        'K10': extract(K10_mat), 'K01': extract(K01_mat),
        'K20': extract(K20_mat), 'K02': extract(K02_mat), 'K11': extract(K11_mat),
        'K30': extract(K30_mat), 'K03': extract(K03_mat),
        'K21': extract(K21_mat), 'K12': extract(K12_mat),
        'Cdd0': extract(Cdd0_mat), 'Cdd1': extract(Cdd1_mat), 'Cdd2': extract(Cdd2_mat),
        'Cgd0': extract(Cgd0_mat)
    }

    print(f"\n=== 本質參數 (VGS={TARGET_VGS*1000}mV, VDS={TARGET_VDS*1000}mV) ===")
    print(f"gm    : {K['K10']*1000:.4f} mS \t gds   : {K['K01']*1e6:.4f} uS")
    print(f"Cdd0  : {K['Cdd0']*1e15:.4f} fF \t Cgd0  : {K['Cgd0']*1e15:.4f} fF")
    print(f"Cdd1  : {K['Cdd1']*1e15:.4f} fF/V (一階電容非線性)")
    print(f"Cdd2  : {K['Cdd2']*1e15:.4f} fF/V^2 (二階電容非線性)")

    # ================= 4. Volterra 級數計算 (複數運算) =================
    w = 2 * np.pi * FIN
    s = 1j * w
    G_L = 1.0 / R_LOAD
    
    # 輸出節點的總線性電容
    C_out0 = K['Cdd0'] + K['Cgd0']
    
    # 建立 Volterra 算子
    # 1. 一階傳輸函數 (包含頻率響應)
    denominator_1 = G_L + K['K01'] + s * C_out0
    a1 = (-K['K10'] + s * K['Cgd0']) / denominator_1
    
    # 2. 二階傳輸函數 (HD2 激勵在 2w)
    # K02 吸收了 Cdd1 帶來的動態位移電流 jw * Cdd1
    K02_HD2 = K['K02'] + 1j * w * K['Cdd1']
    denominator_2 = G_L + K['K01'] + 2 * s * C_out0
    numerator_a2 = K['K20'] + K['K11'] * a1 + K02_HD2 * (a1**2)
    a2 = -numerator_a2 / denominator_2
    
    # 3. 三階傳輸函數 (HD3 激勵在 3w)
    # K03 吸收了 Cdd2 帶來的動態位移電流 jw * Cdd2
    K03_HD3 = K['K03'] + 1j * w * K['Cdd2']
    # 交叉混頻項: w 與 2w 混合，等效係數為 j * 1.5w * Cdd1
    K02_mix = K['K02'] + 1j * 1.5 * w * K['Cdd1']
    
    denominator_3 = G_L + K['K01'] + 3 * s * C_out0
    numerator_a3 = (K['K30'] + K['K21'] * a1 + K['K12'] * (a1**2) + 
                    K03_HD3 * (a1**3) + K['K11'] * a2 + 2 * K02_mix * a1 * a2)
    a3 = -numerator_a3 / denominator_3

    # 計算諧波
    hd2_lin = 0.5 * np.abs(a2 / a1) * V_AMP_IN
    hd3_lin = 0.25 * np.abs(a3 / a1) * (V_AMP_IN**2)
    
    hd2_dbc = 20 * np.log10(hd2_lin)
    hd3_dbc = 20 * np.log10(hd3_lin)
    sfdr_dbc = -max(hd2_dbc, hd3_dbc)
    
    print(f"\n=== 單級 CS 放大器 [高頻 Volterra] 失真 ===")
    print(f"Gain (|a1|) : {np.abs(a1):.4f} V/V (相位: {np.angle(a1, deg=True):.1f}°)")
    print(f"HD2         : {hd2_dbc:.4f} dBc")
    print(f"HD3         : {hd3_dbc:.4f} dBc")
    print(f"SFDR        : {sfdr_dbc:.4f} dBc")

if __name__ == "__main__":
    calculate_volterra_distortion()