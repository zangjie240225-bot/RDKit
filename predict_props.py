#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#预测单个小分子的物理化学性质。使用方法：输入python predict_props.py, 然后回车，提示输入smiles, 输入smiels后回车，即可在屏幕显示这个小分子的理化性质。

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED

def calc_properties(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    props = {
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
    return props


if __name__ == "__main__":
    smiles = input("请输入SMILES字符串: ").strip()
    result = calc_properties(smiles)

    if result:
        print("\n--- 预测结果 ---")
        for k, v in result.items():
            print(f"{k:20s}: {v}")
    else:
        print("❌ 输入的 SMILES 无效，请检查。")

