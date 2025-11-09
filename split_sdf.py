"""
=========================================
脚本名称: split_sdf.py
功能说明: 
    使用 RDKit 将一个大的 SDF 文件拆分成多个小的 SDF 文件。
    每个输出文件包含固定数量的分子（默认 100 个）。

用法:
    python split_sdf.py input.sdf
    python split_sdf.py input.sdf 200
    python split_sdf.py input.sdf 200 output_dir/

示例:
    python split_sdf.py drug_library.sdf 150 split_files/

依赖:
    - Python 3.x
    - RDKit
=========================================
"""

import sys
import os
from rdkit import Chem

# ======= 解析命令行参数 =======
if len(sys.argv) < 2:
    print("用法: python split_sdf.py <输入SDF文件> [每个文件分子数] [输出目录]")
    sys.exit(1)

input_sdf = sys.argv[1]
chunk_size = int(sys.argv[2]) if len(sys.argv) >= 3 else 100
output_dir = sys.argv[3] if len(sys.argv) >= 4 else "."

# ======= 参数检查 =======
if not os.path.exists(input_sdf):
    raise FileNotFoundError(f"❌ 未找到输入文件: {input_sdf}")

os.makedirs(output_dir, exist_ok=True)

# ======= 主程序 =======
mols = Chem.SDMolSupplier(input_sdf)
if not mols:
    raise ValueError("❌ 无法读取分子，请检查 SDF 文件是否有效。")

count = 0
file_index = 1
out_path = os.path.join(output_dir, f"part_{file_index}.sdf")
writer = Chem.SDWriter(out_path)

for mol in mols:
    if mol is None:
        continue
    writer.write(mol)
    count += 1

    if count >= chunk_size:
        writer.close()
        print(f"✅ 已生成 {out_path} ({count} 个分子)")
        file_index += 1
        out_path = os.path.join(output_dir, f"part_{file_index}.sdf")
        writer = Chem.SDWriter(out_path)
        count = 0

writer.close()
print(f"🎯 拆分完成，共生成 {file_index} 个文件。")
