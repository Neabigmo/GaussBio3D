"""
Quick Start Guide for GaussBio3D
GaussBio3D 快速入门指南

This file demonstrates the basic workflow of using GaussBio3D.
本文件演示使用GaussBio3D的基本工作流程。
"""

# 1. Installation / 安装
# =====================
# pip install -e .
# Or install dependencies separately / 或分别安装依赖项:
# pip install numpy biopython rdkit-pypi


# 2. Basic Imports / 基本导入
# ===========================
from gaussbio3d.molecules import Protein, Ligand, NucleicAcid
from gaussbio3d.config import MgliConfig
from gaussbio3d.features.descriptor import global_mgli_descriptor
from gaussbio3d.features.node_features import node_mgli_features
from gaussbio3d.features.pairwise import pairwise_mgli_matrix


# 3. Load Molecules / 加载分子
# =============================

# Load a protein from PDB / 从PDB加载蛋白质
# protein = Protein.from_pdb("path/to/protein.pdb", chain_id="A")

# Load a ligand from SDF / 从SDF加载配体
# ligand = Ligand.from_sdf("path/to/ligand.sdf")

# Load a ligand from SMILES / 从SMILES加载配体
# ligand = Ligand.from_smiles("CC(C)Cc1ccc(cc1)C(C)C(O)=O")

# Load nucleic acid / 加载核酸
# na = NucleicAcid.from_pdb("path/to/dna.pdb", chain_id="D")


# 4. Configure mGLI Parameters / 配置mGLI参数
# ===========================================

# Create configuration / 创建配置
config = MgliConfig(
    # Distance bins (Angstroms) / 距离分箱（埃）
    distance_bins=[0.0, 3.0, 6.0, 10.0, 20.0],
    
    # Use RBF or hard bins / 使用RBF或硬分箱
    use_rbf=False,
    
    # RBF sigma (if use_rbf=True) / RBF sigma（如果use_rbf=True）
    rbf_sigma=None,
    
    # Keep sign or use absolute value / 保留符号或使用绝对值
    signed=False,
    
    # Statistics to compute / 要计算的统计量
    stats=["sum", "mean", "max", "min", "median"],
    
    # Grouping mode for structure A / 结构A的分组模式
    group_mode_A="residue_class",  # or "element" / 或"element"
    
    # Grouping mode for structure B / 结构B的分组模式
    group_mode_B="element",  # or "group" / 或"group"
)


# 5. Compute Features / 计算特征
# ==============================

# Global descriptor (structure pair) / 全局描述符（结构对）
# global_feat = global_mgli_descriptor(protein, ligand, config)
# Shape: (G_A * G_B * K * S,) where:
# 形状：(G_A * G_B * K * S,) 其中：
#   G_A = number of groups in A / A中的组数
#   G_B = number of groups in B / B中的组数
#   K = number of distance scales / 距离尺度数
#   S = number of statistics / 统计量数

# Node-level features / 节点级特征
# prot_node_feat = node_mgli_features(protein, ligand, config)
# Shape: (N_protein, K * S) / 形状：(N_protein, K * S)

# lig_node_feat = node_mgli_features(ligand, protein, config)
# Shape: (N_ligand, K * S) / 形状：(N_ligand, K * S)

# Pairwise matrix / 成对矩阵
# pairwise_mat = pairwise_mgli_matrix(protein, ligand, signed=False)
# Shape: (N_protein, N_ligand) / 形状：(N_protein, N_ligand)


# 6. Task-Specific Helpers / 特定任务辅助函数
# ===========================================

# Drug-Target Interaction (DTI) / 药物-靶点交互
from gaussbio3d.tasks.dti import compute_dti_features
# dti_features = compute_dti_features(
#     pdb_path="protein.pdb",
#     sdf_path="drug.sdf",  # or use smiles="..." / 或使用smiles="..."
#     chain_id="A",
#     config=config
# )
# Returns: {
#   "global_feat": ...,
#   "prot_node_feat": ...,
#   "lig_node_feat": ...,
#   "pairwise_mgli": ...
# }

# Protein-Protein Interaction (PPI) / 蛋白质-蛋白质交互
from gaussbio3d.tasks.ppi import compute_ppi_features
# ppi_features = compute_ppi_features(
#     pdb_path_A="protein_A.pdb",
#     pdb_path_B="protein_B.pdb",
#     chain_id_A="A",
#     chain_id_B="B",
#     config=config
# )

# Molecule-Target Interaction (MTI) / 分子-靶点交互
from gaussbio3d.tasks.mti import compute_mti_features
# mti_features = compute_mti_features(
#     protein_pdb="protein.pdb",
#     na_pdb="dna.pdb",
#     protein_chain="P",
#     na_chain="D",
#     config=config
# )


# 7. Integration with ML Models / 与机器学习模型集成
# =================================================

# Use global features for traditional ML / 将全局特征用于传统机器学习
# from sklearn.ensemble import RandomForestClassifier
# clf = RandomForestClassifier()
# clf.fit(global_features, labels)

# Use node features with GNNs / 将节点特征用于图神经网络
# import torch
# from torch_geometric.nn import GATConv
# # Concatenate mGLI features with other node features
# # 将mGLI特征与其他节点特征连接
# node_feat = torch.cat([plm_embeddings, mgli_features], dim=-1)

# Use pairwise matrix for attention / 将成对矩阵用于注意力机制
# attention_bias = pairwise_mgli_matrix


print("GaussBio3D Quick Start Guide / GaussBio3D 快速入门指南")
print("=" * 60)
print("This file contains code templates for using GaussBio3D.")
print("本文件包含使用GaussBio3D的代码模板。")
print("\nUncomment the code sections to run them with your data.")
print("取消注释代码部分以使用您的数据运行它们。")
print("=" * 60)
