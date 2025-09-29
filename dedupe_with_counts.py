# dedupe_with_counts.py
import pandas as pd
from rdkit import Chem
import sys

if len(sys.argv) < 2:
    print("用法: python dedupe_with_counts.py input.csv")
    sys.exit(1)

input_csv = sys.argv[1]
df = pd.read_csv(input_csv)

# 自动检测 "smiles" 列（忽略大小写）
smiles_col = None
for col in df.columns:
    if col.lower() == "smiles":
        smiles_col = col
        break

if smiles_col is None:
    sys.exit("错误：CSV 里没有找到 SMILES 列，请确认表头")

print(f"检测到分子列: {smiles_col}")

# 生成规范化 SMILES
def canon_smiles(smi):
    mol = Chem.MolFromSmiles(str(smi))
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    return None

df["Canon_SMILES"] = df[smiles_col].apply(canon_smiles)
df = df.dropna(subset=["Canon_SMILES"])

# 统计重复次数
counts = df["Canon_SMILES"].value_counts().reset_index()
counts.columns = ["Canon_SMILES", "Count"]

# 保留第一条作为代表
df_first = df.drop_duplicates(subset=["Canon_SMILES"])
df_dedup = pd.merge(df_first, counts, on="Canon_SMILES")

# 保存
output_csv = input_csv.replace(".csv", "_dedup.csv")
df_dedup.to_csv(output_csv, index=False, encoding="utf-8")

print(f"去重前: {len(df)}, 去重后: {len(df_dedup)}，共去掉 {len(df)-len(df_dedup)} 条重复")
print("前 10 个最常见的分子：")
print(df_dedup.sort_values(by='Count', ascending=False).head(10)[['Canon_SMILES','Count']])
