# ------------------------------------------------------------
# 脚本名称: 计算一个小分子库的diversity
#
# 功能:
#   1. 从 CSV 文件中读取分子 (SMILES 列)
#   2. 生成分子指纹 (Morgan Fingerprint, 类似 ECFP)
#   3. 计算所有分子之间的两两 Tanimoto 相似性
#   4. 输出平均相似性 (mean similarity) 和 多样性 (1 - mean similarity)
#   5. 绘制相似性分布直方图
#   6. 保存 N×N 相似性矩阵到 similarity_matrix.csv
#
# 用法示例:
#   如果 CSV 文件中有一列 "smiles":
#       python calc_diversity.py denovo_results.csv
#
#   如果 SMILES 列的名字不是 "smiles"，比如叫 "SMILES":
#       python calc_diversity.py denovo_results.csv SMILES
#
# 输出:
#   - 终端: 分子数量、平均相似性、多样性
#   - 图表: 相似性分布直方图
#   - 文件: similarity_matrix.csv (相似性矩阵, N×N)
#
# ------------------------------------------------------------

import sys
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, DataStructs
import matplotlib.pyplot as plt

def load_smiles_from_csv(file_path, smiles_col="smiles"):
    """从CSV读取SMILES列表"""
    df = pd.read_csv(file_path)
    if smiles_col not in df.columns:
        raise ValueError(f"❌ CSV文件中没有找到列: {smiles_col}")
    return df[smiles_col].dropna().tolist()

def compute_fingerprints(smiles_list, radius=2, nBits=2048):
    """生成分子指纹"""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fps.append(gen.GetFingerprint(mol))
    return fps

def compute_pairwise_similarities(fps):
    """计算两两相似性"""
    sims = []
    n = len(fps)
    for i in range(n):
        for j in range(i + 1, n):
            sims.append(DataStructs.TanimotoSimilarity(fps[i], fps[j]))
    return sims

def compute_similarity_matrix(fps):
    """生成完整相似性矩阵"""
    n = len(fps)
    sim_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fps[i], fps[j])
    return sim_matrix

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python calc_diversity.py <input_csv> [smiles_col]")
        sys.exit(1)

    input_csv = sys.argv[1]
    smiles_col = sys.argv[2] if len(sys.argv) == 3 else "smiles"

    smiles_list = load_smiles_from_csv(input_csv, smiles_col)
    fps = compute_fingerprints(smiles_list)
    sims = compute_pairwise_similarities(fps)

    mean_sim = np.mean(sims)
    diversity = 1 - mean_sim

    print(f"📊 分子数: {len(smiles_list)}")
    print(f"平均相似性: {mean_sim:.3f}")
    print(f"多样性: {diversity:.3f}")

    # 保存相似性矩阵
    sim_matrix = compute_similarity_matrix(fps)
    output_file = "similarity_matrix.csv"
    pd.DataFrame(sim_matrix, index=smiles_list, columns=smiles_list).to_csv(output_file)
    print(f"✅ 相似性矩阵已保存到 {output_file}")

    # 可视化分布
    plt.hist(sims, bins=20, alpha=0.7, color="blue")
    plt.xlabel("Tanimoto Similarity")
    plt.ylabel("Count")
    plt.title("Pairwise Similarity Distribution")
    plt.show()