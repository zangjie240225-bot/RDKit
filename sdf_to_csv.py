#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
===========================================================
功能：将 SDF 文件转换为 CSV 文件（兼容中文或非 UTF-8 编码）
依赖：RDKit, pandas
-----------------------------------------------------------
用法：
    python sdf_to_csv.py input.sdf output.csv
===========================================================
"""

from rdkit import Chem
import pandas as pd
import sys

input_sdf = sys.argv[1] if len(sys.argv) > 1 else "input.sdf"
output_csv = sys.argv[2] if len(sys.argv) > 2 else "output.csv"

# 关闭 sanitize 以防出错
mols = Chem.SDMolSupplier(input_sdf, sanitize=False, removeHs=False)
data = []

for mol in mols:
    if mol is None:
        continue

    smiles = Chem.MolToSmiles(mol)
    name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
    props = {}

    # 遍历所有属性并安全读取
    for key in mol.GetPropNames():
        try:
            value = mol.GetProp(key)
        except Exception:
            try:
                # 尝试以 GBK 解码原始字节
                value = mol.GetProp(key).encode('latin1').decode('gbk', errors='ignore')
            except Exception:
                value = "?"
        props[key] = value

    props["SMILES"] = smiles
    props["Name"] = name
    data.append(props)

# 保存 CSV，带 BOM 方便 Excel 识别中文
df = pd.DataFrame(data)
df.to_csv(output_csv, index=False, encoding="utf-8-sig")

print(f"✅ 已成功转换 {len(df)} 个分子为 {output_csv}")
