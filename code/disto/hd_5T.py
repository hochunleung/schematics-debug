import numpy as np
import pandas as pd
import re
from scipy.interpolate import RegularGridInterpolator

# ================= 1. 系統配置與 DC 偏壓點 =================
# 請填入 Cadence dcOp 跑出來的實際電壓 (PMOS 請填絕對值)
OP_CONFIG = {
    'M0': {'vgs': 0.3739, 'vds': 0.2216, 'type': 'nch'},      # 尾電流源
    'M2': {'vgs': 0.2284, 'vds': 0.3701, 'type': 'nch_lvt'},  # 輸入對管 (正端)
    'M3': {'vgs': 0.2344, 'vds': 0.2344, 'type': 'nch_lvt'},  # 輸入對管 (負端)
    'M4': {'vgs': 0.3083, 'vds': 0.3083, 'type': 'pch'},      # 負載 (左)
    'M5': {'vgs': 0.3083, 'vds': 0.4440, 'type': 'pch'}       # 負載 (右/輸出)
}

FILE_MAP = {
    'nch_lvt': {'gm': 'nch_lvt_gm.csv', 'gds': 'nch_lvt_gds.csv'},
    'pch':     {'gm': 'pch_gm.csv',     'gds': 'pch_gds.csv'},
    'nch':     {'gm': 'nch_gm.csv',     'gds': 'nch_gds.csv'}
}

# 測試條件
V_AMP_IN = 10e-3  # 1mV
R_LOAD = 1e9     # 假設無顯著外部負載，設為 1M 模擬開路或示波器探棒
# ============================================================

def extract_vgs(col_name):
    match = re.search(r'Vgs[^\d\.\-]*([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name, re.IGNORECASE)
    if match: return float(match.group(1))
    match = re.search(r'([+-]?\d*\.?\d+(?:[eE][+-]?\d+)?)', col_name)
    if match: return float(match.group(1))
    return None

def load_lut(gm_file, gds_file):
    """讀取 CSV 並計算 2D 矩陣的所有偏導項"""
    df_gm = pd.read_csv(gm_file, sep=None, engine='python')
    df_gds = pd.read_csv(gds_file, sep=None, engine='python')
    
    vds_arr = df_gm.iloc[:, 0].astype(float).values
    vgs_arr = np.array([extract_vgs(col) for col in df_gm.columns[1:]])
    
    gm_mat = df_gm.iloc[:, 1:].astype(float).values
    gds_mat = df_gds.iloc[:, 1:].astype(float).values
    
    vds_step = vds_arr[1] - vds_arr[0]
    vgs_step = vgs_arr[1] - vgs_arr[0]
    
    mats = {
        'K10': gm_mat,
        'K01': gds_mat,
        'K20': 0.5 * np.gradient(gm_mat, vgs_step, axis=1),
        'K02': 0.5 * np.gradient(gds_mat, vds_step, axis=0),
        'K11': np.gradient(gm_mat, vds_step, axis=0)
    }
    mats['K30'] = (1/3) * np.gradient(mats['K20'], vgs_step, axis=1)
    mats['K03'] = (1/3) * np.gradient(mats['K02'], vds_step, axis=0)
    mats['K21'] = np.gradient(mats['K20'], vds_step, axis=0)
    mats['K12'] = np.gradient(mats['K02'], vgs_step, axis=1)
    return vds_arr, vgs_arr, mats

def build_Y_matrix(K, RL):
    """建立 3x3 一階線性導納矩陣"""
    Y = np.zeros((3, 3), dtype=float)
    YL = 1.0 / RL
    
    # 節點 0: vs (尾節點)
    Y[0, 0] = K['M0']['K01'] + K['M2']['K10'] + K['M2']['K01'] + K['M3']['K10'] + K['M3']['K01']
    Y[0, 1] = -K['M2']['K01']
    Y[0, 2] = -K['M3']['K10'] - K['M3']['K01']
    
    # 節點 1: vm (鏡像節點)
    Y[1, 0] = -K['M2']['K10'] - K['M2']['K01']
    Y[1, 1] = K['M2']['K01'] + K['M4']['K10'] + K['M4']['K01']
    Y[1, 2] = 0
    
    # 節點 2: vout (輸出節點)
    Y[2, 0] = -K['M3']['K10'] - K['M3']['K01']
    Y[2, 1] = K['M5']['K10']
    Y[2, 2] = K['M3']['K10'] + K['M3']['K01'] + K['M5']['K01'] # + YL
    return Y

def calc_iNL2(K, vg, vd, mos_type):
    """計算進入 Drain 的二階虛擬失真電流"""
    i_nl2 = K['K20']*(vg**2) + K['K11']*(vg*vd) + K['K02']*(vd**2)
    # PMOS 偶次項必須反相
    return -i_nl2 if mos_type == 'pch' else i_nl2

def calc_iNL3(K, vg1, vd1, vg2, vd2, mos_type):
    """計算進入 Drain 的三階虛擬失真電流"""
    i_pure3 = K['K30']*(vg1**3) + K['K21']*(vg1**2)*vd1 + K['K12']*vg1*(vd1**2) + K['K03']*(vd1**3)
    i_mix = 2*K['K20']*vg1*vg2 + K['K11']*(vg1*vd2 + vg2*vd1) + 2*K['K02']*vd1*vd2
    # PMOS：奇次項本質不變，但混頻項是由偶次項衍生，必須反相
    if mos_type == 'pch':
        return i_pure3 - i_mix
    else:
        return i_pure3 + i_mix

def run_volterra_analysis():
    # 1. 自動查表生成 K_dict
    print("=== 正在讀取並萃取 LUT 數據 ===")
    lut_lib = {}
    for tag, files in FILE_MAP.items():
        lut_lib[tag] = load_lut(files['gm'], files['gds'])

    K_dict = {}
    for inst, cfg in OP_CONFIG.items():
        vds_arr, vgs_arr, mats = lut_lib[cfg['type']]
        target = np.array([cfg['vds'], cfg['vgs']])
        K_dict[inst] = {k: RegularGridInterpolator((vds_arr, vgs_arr), m)(target)[0] for k, m in mats.items()}
        print(f"[{inst}] gm: {K_dict[inst]['K10']*1000:.3f} mS | gds: {K_dict[inst]['K01']*1e6:.2f} uS")

    # 2. 求解一階 (線性基頻)
    Y1 = build_Y_matrix(K_dict, R_LOAD)
    gm2 = K_dict['M2']['K10']
    I_in1 = np.array([gm2 * V_AMP_IN, -gm2 * V_AMP_IN, 0])
    V1 = np.linalg.solve(Y1, I_in1)
    
    # 提取各管一階電壓
    v1_s, v1_m, v1_out = V1[0], V1[1], V1[2]
    vg1 = {
        'M0': 0, 'M2': V_AMP_IN - v1_s, 'M3': v1_out - v1_s,
        'M4': v1_m, 'M5': v1_m
    }
    vd1 = {
        'M0': v1_s, 'M2': v1_m - v1_s, 'M3': v1_out - v1_s,
        'M4': v1_m, 'M5': v1_out
    }

    # 3. 求解二階 (二次諧波)
    iNL2 = {inst: calc_iNL2(K_dict[inst], vg1[inst], vd1[inst], OP_CONFIG[inst]['type']) for inst in OP_CONFIG}
    
    I_in2 = np.array([
        iNL2['M2'] + iNL2['M3'] - iNL2['M0'],  # Node vs
        -iNL2['M2'] - iNL2['M4'],              # Node vm
        -iNL2['M3'] - iNL2['M5']               # Node vout
    ])
    
    Y2 = build_Y_matrix(K_dict, R_LOAD) # 靜態下 Y2 與 Y1 相同
    V2 = np.linalg.solve(Y2, I_in2)
    v2_s, v2_m, v2_out = V2[0], V2[1], V2[2]

    # 4. 求解三階 (三次諧波)
    vg2 = {
        'M0': 0, 'M2': 0 - v2_s, 'M3': v2_out - v2_s,
        'M4': v2_m, 'M5': v2_m
    }
    vd2 = {
        'M0': v2_s, 'M2': v2_m - v2_s, 'M3': v2_out - v2_s,
        'M4': v2_m, 'M5': v2_out
    }
    
    iNL3 = {inst: calc_iNL3(K_dict[inst], vg1[inst], vd1[inst], vg2[inst], vd2[inst], OP_CONFIG[inst]['type']) for inst in OP_CONFIG}
    
    I_in3 = np.array([
        iNL3['M2'] + iNL3['M3'] - iNL3['M0'],
        -iNL3['M2'] - iNL3['M4'],
        -iNL3['M3'] - iNL3['M5']
    ])
    
    Y3 = build_Y_matrix(K_dict, R_LOAD)
    V3 = np.linalg.solve(Y3, I_in3)
    
    # 5. 計算失真指標
    a1_closed = V1[2] / V_AMP_IN
#    hd2_lin = np.abs(V2[2] / V1[2])
#    hd3_lin = np.abs(V3[2] / V1[2])

    hd2_lin = 0.5 * np.abs(V2[2] / V1[2])
    hd3_lin = 0.25 * np.abs(V3[2] / V1[2])
    
    print(f"\n=== 5T OTA 閉環靜態失真 (Volterra 矩陣法) ===")
    print(f"閉環增益 (A1): {np.abs(a1_closed):.5f} V/V")
    print(f"HD2        : {20 * np.log10(hd2_lin):.2f} dBc")
    print(f"HD3        : {20 * np.log10(hd3_lin):.2f} dBc")
    print(f"SFDR       : {-max(20 * np.log10(hd2_lin), 20 * np.log10(hd3_lin)):.2f} dBc")

if __name__ == "__main__":
    run_volterra_analysis()