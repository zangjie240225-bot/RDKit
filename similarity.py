
#比较两个小分子的similarity, 使用方法：运行python similarity.py "smiles1" "smiles2"

import sys
import csv
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator

def calc_fingerprint(smiles, radius=2, nBits=2048):
    """新版 Morgan 指纹生成 (推荐)"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)
    return gen.GetFingerprint(mol)

def calc_similarity(smi1, smi2):
    """计算两个分子的相似性（多种度量方法）"""
    fp1 = calc_fingerprint(smi1)
    fp2 = calc_fingerprint(smi2)

    sims = {
        "Tanimoto": DataStructs.TanimotoSimilarity(fp1, fp2),
        "Dice": DataStructs.DiceSimilarity(fp1, fp2),
        "Cosine": DataStructs.CosineSimilarity(fp1, fp2),
        "Sokal": DataStructs.SokalSimilarity(fp1, fp2),
        "Russel": DataStructs.RusselSimilarity(fp1, fp2),
    }
    return sims

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("用法: python similarity.py <SMILES1> <SMILES2>")
        sys.exit(1)

    smi1, smi2 = sys.argv[1], sys.argv[2]
    results = calc_similarity(smi1, smi2)

    # 打印到终端
    print(f"Comparing {smi1} vs {smi2}")
    for k, v in results.items():
        print(f"{k} similarity: {v:.3f}")

    # 保存为 CSV 文件
    output_file = "similarity_results.csv"
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["SMILES1", "SMILES2"] + list(results.keys()))
        writer.writerow([smi1, smi2] + [f"{v:.3f}" for v in results.values()])

    print(f"\n结果已保存到 {output_file}")
