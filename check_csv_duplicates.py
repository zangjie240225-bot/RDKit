# check_csv_duplicates.py
import pandas as pd
from rdkit import Chem
from collections import Counter
import sys

if len(sys.argv) < 2:
    print("用法: python check_csv_duplicates.py input.csv")
    sys.exit(1)

input_csv = sys.argv[1]
df = pd.read_csv(input_csv)

# 自动识别 SMILES 列（忽略大小写）
smiles_col = None
for col in df.columns:
    if col.lower() == "smiles":
        smiles_col = col
        break

if smiles_col is None:
    sys.exit("错误：CSV 文件中没有找到 SMILES 列")

print(f"检测到分子列: {smiles_col}")

# 生成规范化 SMILES
def canon_smiles(smi):
    mol = Chem.MolFromSmiles(str(smi))
    if mol:
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    return None

df["Canon_SMILES"] = df[smiles_col].apply(canon_smiles)
df = df.dropna(subset=["Canon_SMILES"])

# 统计
total = len(df)
counts = Counter(df["Canon_SMILES"])
unique = len(counts)
dups = total - unique

print(f"总分子数: {total}")
print(f"唯一分子数: {unique}")
print(f"重复数: {dups}")

if dups > 0:
    print("\n前 10 个最常见的重复分子：")
    for smi, c in counts.most_common(10):
        if c > 1:
            print(f"{smi}  → {c} 次")
else:
    print("\n没有检测到重复分子。")
