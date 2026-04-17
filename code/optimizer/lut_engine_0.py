import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator, interp1d

class MosData:
    def __init__(self, csv_path, W_ref=4e-6, is_pmos=False):
        """
        加載並處理 MOSFET LUT 數據。
        :param csv_path: CSV 檔案路徑
        :param W_ref: 提取數據時使用的單根 Finger 寬度 (預設 2*2um)
        :param is_pmos: 是否為 PMOS (自動對電壓與電流取絕對值)
        """
        print(f"Loading {csv_path}...")
        self.df = pd.read_csv(csv_path)
        
        # 1. 處理 PMOS 的負值
        if is_pmos:
            self.df = self.df.abs()
            
        # 2. 過濾極端低電流 (避免除以零)，設定一個微小的電流下限
        epsilon = 1e-15
        self.df['ID'] = self.df['ID'].clip(lower=epsilon)
        self.df['GDS'] = self.df['GDS'].clip(lower=epsilon)
        self.df['GM'] = self.df['GM'].clip(lower=epsilon)
        
        # 3. 計算密度與本質參數
        self.df['ID_W'] = self.df['ID'] / W_ref
        self.df['GM_ID'] = self.df['GM'] / self.df['ID']
        self.df['GM_GDS'] = self.df['GM'] / self.df['GDS'] # 本質增益 (Self-Gain)
        self.df['CGS_W'] = self.df['CGS'] / W_ref
        self.df['CGD_W'] = self.df['CGD'] / W_ref
        self.df['CDB_W'] = self.df['CDB'] / W_ref
        self.df['CGG_W'] = self.df['CGG'] / W_ref
        
        # 估算轉移頻率 fT (大約等於 gm / (2 * pi * Cgg))
        self.df['fT'] = self.df['GM'] / (2 * np.pi * self.df['CGG'].clip(lower=1e-18))
        
        # 4. 提取唯一的座標軸 (保證排序)
        self.L_axis = np.sort(self.df['L'].unique())
        self.VGS_axis = np.sort(self.df['VGS'].unique())
        self.VDS_axis = np.sort(self.df['VDS'].unique())
        
        # 5. 建立 3D 網格字典 (以 L, VGS, VDS 為軸)
        self.grids = {}
        # 將 DataFrame 設置多重索引，並轉換為 3D Numpy Array
        df_indexed = self.df.set_index(['L', 'VGS', 'VDS']).sort_index()
        
        # 檢查網格是否完整
        expected_size = len(self.L_axis) * len(self.VGS_axis) * len(self.VDS_axis)
        if len(df_indexed) != expected_size:
            print(f"Warning: Grid size mismatch! Expected {expected_size}, got {len(df_indexed)}. Interpolation may fail.")

        # [階段一] 先載入基礎參數並 reshape 成 3D 陣列
        base_params = ['ID_W', 'GM_ID', 'GM_GDS', 'fT', 'GM', 'GDS', 'CGS_W', 'CGD_W', 'CDB_W', 'CGG_W']
        for param in base_params:
            self.grids[param] = df_indexed[param].values.reshape(
                (len(self.L_axis), len(self.VGS_axis), len(self.VDS_axis))
            )

        # [階段二] 提取 3D 陣列，計算非線性高階導數 Kij
        vgs_step = self.VGS_axis[1] - self.VGS_axis[0]
        vds_step = self.VDS_axis[1] - self.VDS_axis[0]
        # print(f"vgs_step: {vgs_step}")
        # print(f"vds_step: {vgs_step}")
        
        # 將 GM 和 GDS 標準化為單位寬度 (W_ref) 的導數
        gm_w = self.grids['GM'] / W_ref
        gds_w = self.grids['GDS'] / W_ref
        
        self.grids['K10_W'] = gm_w
        self.grids['K01_W'] = gds_w
        self.grids['K20_W'] = 0.5 * np.gradient(gm_w, vgs_step, axis=1)
        self.grids['K02_W'] = 0.5 * np.gradient(gds_w, vds_step, axis=2)
        self.grids['K11_W'] = np.gradient(gm_w, vds_step, axis=2)
        self.grids['K30_W'] = (1/3) * np.gradient(self.grids['K20_W'], vgs_step, axis=1)
        self.grids['K03_W'] = (1/3) * np.gradient(self.grids['K02_W'], vds_step, axis=2)
        self.grids['K21_W'] = np.gradient(self.grids['K20_W'], vds_step, axis=2)
        self.grids['K12_W'] = np.gradient(self.grids['K02_W'], vgs_step, axis=1)

        # 整合所有需要建立內插器的參數名單
        all_params = base_params + ['K10_W', 'K01_W', 'K20_W', 'K02_W', 'K11_W', 'K30_W', 'K03_W', 'K21_W', 'K12_W']

        # [階段三] 建立 3D 內插器 (For Forward Lookup)
        self.interps = {}
        for param in all_params:
            self.interps[param] = RegularGridInterpolator(
                (self.L_axis, self.VGS_axis, self.VDS_axis), 
                self.grids[param], 
                bounds_error=False, fill_value=None # 允許輕微外推
            )

    def lookup_vgs_by_idw(self, target_idw, L, VDS_guess):
        """
        給定目標電流密度 (ID/W), L, VDS, 反推需要的 VGS。
        """
        vgs_array = self.VGS_axis
        # 構造 3D 座標點陣列
        pts = np.column_stack((np.full_like(vgs_array, L), vgs_array, np.full_like(vgs_array, VDS_guess)))
        
        # 取出該 L 與 VDS 下的 ID_W 曲線
        idw_curve = self.interps['ID_W'](pts)
        
        # 過濾有效單調遞增區間
        valid_idx = idw_curve > 1e-12
        idw_curve_valid = idw_curve[valid_idx]
        vgs_valid = vgs_array[valid_idx]
        
        # 建立 1D 插值器 (X: ID_W, Y: VGS)
        inv_interp = interp1d(idw_curve_valid, vgs_valid, kind='linear', bounds_error=False, fill_value="extrapolate")
        
        return float(inv_interp(target_idw))

    def lookup_by_vgs(self, param, L, VGS, VDS):
        """正向查表：給定 L, VGS, VDS 查參數"""
        pts = np.array([[L, VGS, VDS]])
        return self.interps[param](pts)[0]

    def lookup_by_gmid(self, target_gmid, L, VDS):
        """
        逆向查表 (優化器核心)：給定目標 gm/ID, L, VDS，反推 VGS 與所有參數。
        這使用了降維插值法。
        """
        # 1. 在指定的 L 和 VDS 下，抽出該 1D 切片的 VGS 和 GM_ID 曲線
        # 我們直接對 VGS 軸上所有的點做一次插值
        vgs_array = self.VGS_axis
        pts = np.column_stack((np.full_like(vgs_array, L), vgs_array, np.full_like(vgs_array, VDS)))
        
        gmid_curve = self.interps['GM_ID'](pts)
        
        # 2. 過濾掉 GM_ID 曲線中非單調遞減的部分 (通常在 VGS 極小且 ID 接近雜訊時發生)
        # 確保 gmid_curve 是嚴格遞減的，這樣才能安全地反向插值
        valid_idx = gmid_curve > 0
        gmid_curve = gmid_curve[valid_idx]
        vgs_valid = vgs_array[valid_idx]
        
        # 將曲線反轉 (從嚴格遞減變成嚴格遞增)，讓 interp1d 可以工作
        gmid_curve_rev = gmid_curve[::-1]
        vgs_valid_rev = vgs_valid[::-1]
        
        # 3. 建立 1D 逆向插值器：輸入 gm/ID，輸出 VGS
        inv_interp = interp1d(gmid_curve_rev, vgs_valid_rev, kind='linear', bounds_error=False, fill_value="extrapolate")
        
        vgs_target = float(inv_interp(target_gmid))
        
        # 4. 有了 VGS，就可以正向查出所有的參數密度了！
        result = {'VGS': vgs_target}
        for param in ['ID_W', 'GM_GDS', 'CGS_W', 'CGD_W', 'CDB_W']:
            result[param] = self.lookup_by_vgs(param, L, vgs_target, VDS)
            
        return result

# =====================================================================
# 可視化檢查代碼 (Data Verification)
# =====================================================================
def plot_mos_characteristics(mos_data, title="MOSFET Characteristics"):
    plt.figure(figsize=(12, 8))
    
    # 選擇要畫的 L 列表 (取頭、中、尾)
    L_plot = [mos_data.L_axis[0], mos_data.L_axis[len(mos_data.L_axis)//2], mos_data.L_axis[-1]]
    VDS_plot = mos_data.VDS_axis[-1] # 看飽和區特性
    
    # 1. gm/ID vs VGS
    plt.subplot(2, 2, 1)
    for L in L_plot:
        vgs = mos_data.VGS_axis
        gmid = mos_data.interps['GM_ID'](np.column_stack((np.full_like(vgs, L), vgs, np.full_like(vgs, VDS_plot))))
        plt.plot(vgs, gmid, label=f'L={L*1e9:.0f}nm')
    plt.xlabel('VGS (V)'); plt.ylabel('gm/ID (V^-1)'); plt.title('Efficiency vs Bias')
    plt.grid(True); plt.legend()

    # 2. ID/W vs VGS (Log scale)
    plt.subplot(2, 2, 2)
    for L in L_plot:
        vgs = mos_data.VGS_axis
        idw = mos_data.interps['ID_W'](np.column_stack((np.full_like(vgs, L), vgs, np.full_like(vgs, VDS_plot))))
        plt.semilogy(vgs, idw, label=f'L={L*1e9:.0f}nm')
    plt.xlabel('VGS (V)'); plt.ylabel('ID/W (A/m)'); plt.title('Subthreshold Behavior')
    plt.grid(True); plt.legend()

    # 3. gm/gds (Intrinsic Gain) vs gm/ID
    plt.subplot(2, 2, 3)
    for L in L_plot:
        vgs = mos_data.VGS_axis
        pts = np.column_stack((np.full_like(vgs, L), vgs, np.full_like(vgs, VDS_plot)))
        gmid = mos_data.interps['GM_ID'](pts)
        gain = mos_data.interps['GM_GDS'](pts)
        plt.plot(gmid, 20*np.log10(gain), label=f'L={L*1e9:.0f}nm')
    plt.xlabel('gm/ID (V^-1)'); plt.ylabel('Intrinsic Gain (dB)'); plt.title('Gain vs Efficiency')
    plt.xlim(5, 25); plt.grid(True); plt.legend()

    # 4. fT vs gm/ID
    plt.subplot(2, 2, 4)
    for L in L_plot:
        vgs = mos_data.VGS_axis
        pts = np.column_stack((np.full_like(vgs, L), vgs, np.full_like(vgs, VDS_plot)))
        gmid = mos_data.interps['GM_ID'](pts)
        ft = mos_data.interps['fT'](pts)
        plt.plot(gmid, ft/1e9, label=f'L={L*1e9:.0f}nm')
    plt.xlabel('gm/ID (V^-1)'); plt.ylabel('fT (GHz)'); plt.title('Transit Frequency vs Efficiency')
    plt.xlim(5, 25); plt.grid(True); plt.legend()

    plt.suptitle(f"{title} @ VDS={VDS_plot:.2f}V", fontsize=16)
    plt.tight_layout()
    plt.show()

# 執行測試
if __name__ == "__main__":
    # 替換為你的路徑
    nch_lvt = MosData('nch_lvt_lut.csv', W_ref=4e-6, is_pmos=False)
    
    # 可視化檢查數據品質
    plot_mos_characteristics(nch_lvt, title="NCH_LVT LUT Verification")
    
    # 測試優化器會用到的「逆向查表」函數
    print("\n--- Optimizer Lookup Test ---")
    res = nch_lvt.lookup_by_gmid(target_gmid=15.0, L=150e-9, VDS=0.9)
    print(f"Target: gm/ID = 15.0, L = 150nm, VDS = 0.9V")
    print(f"-> Required VGS : {res['VGS']:.4f} V")
    print(f"-> ID/W         : {res['ID_W']:.2f} uA/um")
    print(f"-> Intrinsic Gain: {20*np.log10(res['GM_GDS']):.1f} dB")