#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
把分子生成程序的 CSV 转成 SDF（2D、隐式氢，适合后续在 Maestro 里做 LigPrep）
改进点：确保输出不包含显式氢（双重移除与检查）
用法: python convert_csv_to_sdf_no_explicitH.py input.csv
"""

import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdDepictor

# 参数检查
if len(sys.argv) != 2:
    print("用法: python convert_csv_to_sdf_no_explicitH.py input.csv")
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
warn_count = 0

for i, smi in enumerate(df[smiles_col]):
    if pd.isna(smi):
        continue
    smi = str(smi).strip()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"⚠️ 第{i}行 SMILES 无效: {smi}")
        continue

    # 先更新属性缓存并标准化
    try:
        mol.UpdatePropertyCache(strict=False)
    except Exception:
        pass

    # 彻底移除任何显式氢（双重保险）
    try:
        # RemoveHs 返回新的分子引用或就地修改（视版本），我们重新赋值以确保使用最新对象
        mol = Chem.RemoveHs(mol)
    except Exception:
        # 如果 RemoveHs 失败，尝试遍历删 H（极少出现）
        mol = Chem.RWMol(mol)
        to_remove = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'H']
        for idx in sorted(to_remove, reverse=True):
            mol.RemoveAtom(idx)
        mol = mol.GetMol()

    # 再次 sanitize，防止奇怪的 valence/explicit H 状态
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # 若 sanitize 失败，尝试用较宽松参数
        try:
            Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
        except Exception:
            pass

    # 检查是否仍有显式氢（用于 debug/提示）
    nH = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'H')
    if nH > 0:
        warn_count += 1
        print(f"⚠️ 警告：第{i}行仍检测到 {nH} 个显式氢，将再次移除。SMILES: {smi}")
        # 再次尝试移除
        mol = Chem.RemoveHs(mol)
        # 尝试 sanitize 再次确认
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass

    # 生成 2D 坐标（不加氢，保持隐式氢）
    try:
        rdDepictor.Compute2DCoords(mol)
    except Exception as e:
        print(f"⚠️ 第{i}行生成2D坐标失败: {smi}，错误: {e}")
        continue

    # 名称（优先 name/id/title 列）
    name = None
    for col in df.columns:
        if col.lower() in ["name", "id", "title"]:
            val = df.at[i, col]
            if not pd.isna(val):
                name = str(val)
                break
    mol.SetProp("_Name", name if name else f"mol_{i}")

    # 其他列写入 SDF 属性（不写 smiles 列）
    for col in df.columns:
        if col == smiles_col:
            continue
        val = df.at[i, col]
        if pd.isna(val):
            continue
        # 属性值必须是字符串
        mol.SetProp(str(col), str(val))

    # 最终检查：确保写入前分子没有显式氢
    final_nH = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == 'H')
    if final_nH > 0:
        print(f"❗ 最终仍检测到显式氢（{final_nH}）在第{i}行，SMILES: {smi}，该分子将被跳过以避免显式氢写入SDF。")
        continue

    writer.write(mol)
    count += 1

writer.close()
print(f"✅ 成功导出 {count} 个分子到 {output_sdf}（2D、隐式氢，建议后续用 LigPrep 生成3D并加氢）")
if warn_count:
    print(f"提示：有 {warn_count} 个分子在处理中曾出现显式氢并被移除。")
