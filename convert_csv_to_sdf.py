#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把分子生成程序的 CSV 转成 SDF（2D、隐式氢，适合后续在 Maestro 里做 LigPrep）
用法: python convert_csv_to_sdf.py input.csv。缺点是只能识别一列的smiles。
"""

import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDepictor

# 参数检查
if len(sys.argv) != 2:
    print("用法: python convert_csv_to_sdf.py input.csv")
    sys.exit(1)

input_csv = sys.argv[1]
if not os.path.exists(input_csv):
    print(f"❌ 找不到输入文件: {input_csv}")
    sys.exit(1)

# 输出文件名同名 .sdf
output_sdf = os.path.splitext(input_csv)[0] + ".sdf"

# 读 CSV
df = pd.read_csv(input_csv)

# 自动检测 SMILES 列（支持 'smiles' 或 'generated_smiles' 等大小写变体）
smiles_col = None
for col in df.columns:
    if "smiles" in col.lower():
        smiles_col = col
        break

if smiles_col is None:
    print("❌ 没找到 SMILES 列，请确认 CSV 里含有 'smiles' 相关列名")
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

    # 生成 2D 坐标（不加氢，保持隐式氢；更适合在 Maestro 里 LigPrep）
    try:
        rdDepictor.Compute2DCoords(mol)
    except Exception as e:
        print(f"⚠️ 第{i}行生成2D坐标失败: {smi}，错误: {e}")
        continue

    # 可选：如 SMILES 中包含显式氢（例如 [H]N），去掉它们以避免 ChemDraw/导入全显氢
    try:
        mol = Chem.RemoveHs(mol)
    except Exception:
        pass

    # 名称
    name = None
    for col in df.columns:
        if col.lower() in ["name", "id", "title"]:
            name = str(df.at[i, col])
            break
    mol.SetProp("_Name", name if name else f"mol_{i}")

    # 其他列写入 SDF 属性
    for col in df.columns:
        if col == smiles_col:
            continue
        val = df.at[i, col]
        if pd.isna(val):
            continue
        mol.SetProp(col, str(val))

    writer.write(mol)
    count += 1

writer.close()
print(f"✅ 成功导出 {count} 个分子到 {output_sdf}（2D、隐式氢，建议后续用 LigPrep 生成3D并加氢）")




