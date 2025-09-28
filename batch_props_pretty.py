#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#批量预测小分子的性质，输入文件为smi格式的文本文档 （每行一个smile)。使用方法：输入python batch_props_pretty.py input.smi output.csv (其中input.smi为输入的smi文本文件，output.csv为输出文件)
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED


def calc_properties(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    return {
        "SMILES": smiles,
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
    """数值保留 2 位小数，字符串保持不变"""
    if isinstance(val, (int, float)):
        return round(val, 2)
    return val


def main(input_file, output_file):
    data = []
    with open(input_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            smiles = parts[0]
            name = parts[1] if len(parts) > 1 else ""
            props = calc_properties(smiles)
            if props:
                props["Name"] = name
                # 格式化数值
                props = {k: format_value(v) for k, v in props.items()}
                data.append(props)

    df = pd.DataFrame(data)
    # 调整列顺序（名字放最前面）
    cols = ["Name", "SMILES"] + [c for c in df.columns if c not in ["Name", "SMILES"]]
    df = df[cols]

    df.to_csv(output_file, index=False, encoding="utf-8-sig")
    print(f"✅ 已保存到 {output_file}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python batch_props_pretty.py input.smi output.csv")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])

