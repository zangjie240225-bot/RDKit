多样性计算
similarity: 计算两个小分子的相似性，比如Tanimoto参数
calc_diversity: 计算一个小分子库的多样性

理化性质计算
predict_props: 计算一个小分子的各项理化性质，输入为小分子的smiles
batch_props_pretty: 批量计算小分子的理化性质，输入为.smi文件
batch_props_from_csv.py: 批量计算csv文件中小分子的理化性质，输入为.csv文件。需要定义输入和输出文件

库条件筛选
filter_library: 通过RDKit本身的计算，只保留一个小分子库中，满足制定条件的分子。输入为.csv文件
filter_wOwnProp.py:有些下载下来的文件，本身就包含了计算好的理化性质，使用这些自带的理化性质，进行满足定制条件的筛选

格式转换
convert_csv_to_sdf: 把csv文件转换成mastro对接可用的sdf文件
sdf_to_csv: 把sdf文件转换成csv文件，有时候sdf文件用maestro打开会出现错误。
库查重，去重
check_csv_duplicates: 检查一个库文件csv中的重复分子数，基于AI的分子生成会产生这种情况
dedupe_with_counts:对库文件csv进行去重，并统计前10个最大重复的

sdf文件拆分
split_sdf.py：将一个sdf文件按照分子的数量拆分成很多份。可以在命令中定义每份包含多少个分子
