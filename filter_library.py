# ------------------------------------------------------------
# 脚本名称: filter_library.py
#
# 功能:
#   1. 从 CSV 文件中读取分子 (自动识别 SMILES 列)
#   2. 自动兼容 Excel 导出的 CSV (首行可能是 "sep=,")
#   3. 用 RDKit 计算理化性质 (MolWt, LogP, TPSA, HBD, HBA)
#   4. 按条件筛选分子 (所有参数都有最小值和最大值)
#   5. 输出新的 CSV 文件 (filtered_results.csv)
#
# 用法示例:
#   python filter_library.py input.csv
#
# ------------------------------------------------------------

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

# 筛选条件 (可修改范围)
MW_MIN, MW_MAX = 240, 300       # 分子量范围
LOGP_MIN, LOGP_MAX = 2, 3       # cLogP 范围
TPSA_MIN, TPSA_MAX = 0, 70      # 极性表面积
HBD_MIN, HBD_MAX = 1, 2         # 氢键供体数
HBA_MIN, HBA_MAX = 0, 10        # 氢键受体数

def calc_props(smiles):
    """计算分子性质"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return {
        "MolWt": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": Descriptors.TPSA(mol),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
    }

def pass_filters(props):
    """检查分子是否符合筛选条件"""
    if not props:
        return False
    if not (MW_MIN <= props["MolWt"] <= MW_MAX):
        return False
    if not (LOGP_MIN <= props["LogP"] <= LOGP_MAX):
        return False
    if not (TPSA_MIN <= props["TPSA"] <= TPSA_MAX):
        return False
    if not (HBD_MIN <= props["HBD"] <= HBD_MAX):
        return False
    if not (HBA_MIN <= props["HBA"] <= HBA_MAX):
        return False
    return True

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python filter_library.py <input_csv>")
        sys.exit(1)

    input_csv = sys.argv[1]

    # 自动检测是否有 "sep=," 行
    with open(input_csv) as f:
        first_line = f.readline().strip()
    skiprows = 1 if first_line.lower().startswith("sep=") else 0

    # 读取 CSV
    df = pd.read_csv(input_csv, skiprows=skiprows)

    # 自动检测 SMILES 列（允许大小写/变体）
    possible_cols = [c for c in df.columns if "SMI" in c.upper()]
    if not possible_cols:
        raise ValueError(f"❌ 输入文件中没有找到 SMILES 列，请检查表头: {list(df.columns)}")
    smiles_col = possible_cols[0]
    print(f"✅ 自动检测到 SMILES 列: {smiles_col}")

    results = []
    for smi in df[smiles_col].dropna():
        props = calc_props(smi)
        if pass_filters(props):
            results.append({"smiles": smi, **props})

    df_out = pd.DataFrame(results)
    output_file = "filtered_results.csv"
    df_out.to_csv(output_file, index=False)

    print(f"✅ 输入分子数: {len(df)}")
    print(f"✅ 通过筛选分子数: {len(df_out)}")
    print(f"📂 结果已保存到 {output_file}")