import pandas as pd
import numpy as np

print("Loading CSV files...")
df_dc = pd.read_csv("nch_lvt_dc_lut.csv")
df_noise = pd.read_csv("../noise/nch_lvt_noise_lut.csv")

def create_robust_keys(df):
    df['L_key'] = (df['L'] * 1e12).round().astype(int)   # 轉為皮米 (pm)
    df['VGS_key'] = (df['VGS'] * 1e6).round().astype(int) # 轉為微伏 (uV)
    df['VDS_key'] = (df['VDS'] * 1e6).round().astype(int) # 轉為微伏 (uV)
    return df

print("Preprocessing keys...")
df_dc = create_robust_keys(df_dc)
df_noise = create_robust_keys(df_noise)

# 執行合併 (使用我們建立的整數 Key)
# 1. 執行合併前，把右表 (df_noise) 裡會重複的原始欄位先刪除
# 這樣 Pandas 就不會產生 _x 和 _y 的後綴了
cols_to_drop = ['L', 'VGS', 'VDS']
df_noise_clean = df_noise.drop(columns=cols_to_drop)

# 2. 執行合併 (使用我們建立的整數 Key)
print("Merging dataframes...")
merge_keys = ['L_key', 'VGS_key', 'VDS_key']
df_merged = pd.merge(df_dc, df_noise_clean, on=merge_keys, how='inner')

# 3. 移除臨時建立的整數 Key
df_merged = df_merged.drop(columns=merge_keys)
# df_merged = df_merged.loc[:,~df_merged.columns.duplicated()]
# (後續的驗證行數和 to_csv 保持不變)

# 4. 驗證行數 (非常重要！)
expected_count = 816310
actual_count = len(df_merged)
print(f"Merged count: {actual_count}")

if actual_count == expected_count:
    print("Success! Row count matches perfectly.")
else:
    print(f"Warning! Row count mismatch. Expected {expected_count}, got {actual_count}.")
    # 如果還是不對，可能是原始資料就有重複，強制去重
    df_merged = df_merged.drop_duplicates(subset=['L', 'VGS', 'VDS'])
    print(f"Count after forced deduplication: {len(df_merged)}")

# 5. 輸出
df_merged.to_csv("nch_lvt_lut.csv", index=False, float_format='%.6e')
print("Saved to nch_lvt_lut.csv")