#把分子生成程序生成的分子文件xx.csv文件中smiles那一行的分子转换成sdf文件。 使用方法：python convert_csv_to_sdf.py input.csv 


import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

if len(sys.argv) != 2:
    print("用法: python convert_csv_to_sdf.py input.csv ")
    sys.exit(1)

input_csv = sys.argv[1]
if not os.path.exists(input_csv):
    print(f"❌ 找不到输入文件: {input_csv}")
    sys.exit(1)

# 自动生成输出文件名
output_sdf = os.path.splitext(input_csv)[0] + ".sdf"

# 读CSV
df = pd.read_csv(input_csv)

# 自动检测 SMILES 列
smiles_col = None
for col in df.columns:
    if "smiles" in col.lower():
        smiles_col = col
        break

if smiles_col is None:
    print("❌ 没找到 SMILES 列，请确认 CSV 里有 'smiles' 列")
    sys.exit(1)

writer = Chem.SDWriter(output_sdf)
count = 0

for i, smi in enumerate(df[smiles_col]):
    if pd.isna(smi):
        continue
    smi = str(smi).strip()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        continue
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    # 如果有名字列就用名字，否则用编号
    name = None
    for col in df.columns:
        if col.lower() in ["name", "id", "title"]:
            name = str(df[col].iloc[i])
            break
    mol.SetProp("_Name", name if name else f"mol_{i}")

    writer.write(mol)
    count += 1

writer.close()
print(f"✅ 成功导出 {count} 个分子到 {output_sdf}")

