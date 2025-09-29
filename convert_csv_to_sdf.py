# 把分子生成程序生成的分子文件 xx.csv 文件中 smiles 那一列的分子转换成 sdf 文件。
# 使用方法：python convert_csv_to_sdf.py input.csv
# 不需要写输出文件名称和格式。

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
        print(f"⚠️ 第{i}行 SMILES 无效: {smi}")
        continue

    mol = Chem.AddHs(mol)

    # 尝试生成构象
    res = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    if res != 0:
        print(f"⚠️ 第{i}行无法生成3D构象: {smi}")
        continue

    try:
        AllChem.UFFOptimizeMolecule(mol)
    except Exception as e:
        print(f"⚠️ 第{i}行优化失败: {smi}, 错误: {e}")
        continue

    # 如果有名字列就用名字，否则用编号
    name = None
    for col in df.columns:
        if col.lower() in ["name", "id", "title"]:
            name = str(df[col].iloc[i])
            break
    mol.SetProp("_Name", name if name else f"mol_{i}")

    # 把 CSV 里其他列写入 SDF 属性
    for col in df.columns:
        if col != smiles_col:
            mol.SetProp(col, str(df[col].iloc[i]))

    writer.write(mol)
    count += 1

writer.close()
print(f"✅ 成功导出 {count} 个分子到 {output_sdf}")



