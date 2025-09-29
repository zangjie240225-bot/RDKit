# ------------------------------------------------------------
# 脚本名称: filter_library_existing.py
#
# 功能:
#   1. 从 CSV 文件中读取分子 (自动兼容 Excel 导出，首行 "sep=,")
#   2. 使用 CSV 文件中已有的理化性质列 (MW, ClogP, TPSA, HBD, HBA, logS, RotBonds)
#   3. 按条件筛选分子，保留原始所有列 (除非用户另行删除)
#   4. 输出新的 CSV 文件 (filtered_results.csv)，并保留首行 "sep=,"
#
# 用法示例:
#   python filter_library_existing.py input.csv
#
# ------------------------------------------------------------

import sys
import pandas as pd

# ---------------- 筛选条件 (可修改范围) ----------------
MW_MIN, MW_MAX = 240, 300                 # 分子量范围
MW_DESALTED_MIN, MW_DESALTED_MAX = 240, 300  # 脱盐分子量范围
LOGP_MIN, LOGP_MAX = 2, 3                 # cLogP 范围
LOGS_MIN, LOGS_MAX = -10, 10              # logS 范围
TPSA_MIN, TPSA_MAX = 0, 70                # 极性表面积范围
HBD_MIN, HBD_MAX = 1, 2                   # 氢键供体数范围
HBA_MIN, HBA_MAX = 0, 10                  # 氢键受体数范围
ROTBONDS_MIN, ROTBONDS_MAX = 0, 10        # 可旋转键数范围
# ------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python filter_library_existing.py <input_csv>")
        sys.exit(1)

    input_csv = sys.argv[1]

    # 检查第一行是否 "sep=,"
    with open(input_csv) as f:
        first_line = f.readline().strip()
    skiprows = 1 if first_line.lower().startswith("sep=") else 0

    # 读取 CSV（跳过 sep=,）
    df = pd.read_csv(input_csv, skiprows=skiprows)

    # 优先选择 MW (desalted)，否则用 MW
    mw_col = "MW (desalted)" if "MW (desalted)" in df.columns else "MW"

    # 筛选条件
    df_filtered = df[
        (df["MW"] >= MW_MIN) & (df["MW"] <= MW_MAX) &
        (df[mw_col] >= MW_DESALTED_MIN) & (df[mw_col] <= MW_DESALTED_MAX) &
        (df["ClogP"] >= LOGP_MIN) & (df["ClogP"] <= LOGP_MAX) &
        (df["logS"] >= LOGS_MIN) & (df["logS"] <= LOGS_MAX) &
        (df["TPSA"] >= TPSA_MIN) & (df["TPSA"] <= TPSA_MAX) &
        (df["HBD"] >= HBD_MIN) & (df["HBD"] <= HBD_MAX) &
        (df["HBA"] >= HBA_MIN) & (df["HBA"] <= HBA_MAX) &
        (df["RotBonds"] >= ROTBONDS_MIN) & (df["RotBonds"] <= ROTBONDS_MAX)
    ]

    # 输出并保留 sep=,
    output_file = "filtered_results.csv"
    with open(output_file, "w") as f:
        f.write("sep=,\n")  # 保留第一行
    df_filtered.to_csv(output_file, mode="a", index=False)

    print(f"✅ 输入分子数: {len(df)}")
    print(f"✅ 通过筛选分子数: {len(df_filtered)}")
    print(f"📂 结果已保存到 {output_file}")