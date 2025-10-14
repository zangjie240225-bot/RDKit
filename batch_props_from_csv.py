#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量预测小分子的理化性质
输入文件为 CSV 格式，其中必须包含一列名为 "SMILES"
使用方法:
    python batch_props_from_csv.py input.csv output.csv
"""

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED


def calc_properties(smiles: str):
    """计算分子性质"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "分子量 (MolWt)": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "氢键供体数 (HBD)": Lipinski.NumHDonors(mol),
        "氢键受体数 (HBA)": Lipinski.NumHAcceptors(mol),
        "极性表面积 (TPSA)": Descriptors.TPSA(mol),
        "可旋转键数": Lipinski.NumRotatableBonds(mol),
        "芳香环数量": Descriptors.NumAromaticRings(mol),
        "重原子数": Descriptors.HeavyAtomCount(mol),
        "QED (成药性评分)": QED.qed(mol),
    }


def format_value(val):
    """数值保留 2 位小数，字符串保持原样"""
    if isinstance(val, (int, float)):
        return round(val, 2)
    return val


def main(input_csv, output_csv):
    df = pd.read_csv(input_csv)
    if "SMILES" not in df.columns:
        raise ValueError("输入文件必须包含一列名为 'SMILES'")

    results = []
    for i, row in df.iterrows():
        smiles = str(row["SMILES"]).strip()
        if not smiles:
            continue

        props = calc_properties(smiles)
        if props:
            # 合并原始信息 + 性质
            props.update(row.to_dict())
            props = {k: format_value(v) for k, v in props.items()}
            results.append(props)

    df_out = pd.DataFrame(results)

    # 确保 SMILES 和原始列在前面
    base_cols = [c for c in df.columns if c in df_out.columns]
    prop_cols = [c for c in df_out.columns if c not in base_cols]
    df_out = df_out[base_cols + prop_cols]

    df_out.to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"✅ 已保存到 {output_csv}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python batch_props_from_csv.py input.csv output.csv")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])
