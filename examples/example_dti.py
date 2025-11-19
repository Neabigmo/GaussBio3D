"""
Example: Basic DTI Feature Computation
示例：基本DTI特征计算

This script demonstrates how to use GaussBio3D to compute mGLI features
for a drug-target interaction pair.

本脚本演示如何使用GaussBio3D计算药物-靶点交互对的mGLI特征。
"""

import numpy as np
from gaussbio3d.molecules import Protein, Ligand
from gaussbio3d.config import MgliConfig
from gaussbio3d.features.descriptor import global_mgli_descriptor
from gaussbio3d.features.node_features import node_mgli_features
from gaussbio3d.features.pairwise import pairwise_mgli_matrix
from gaussbio3d.tasks.dti import compute_dti_features


def example_basic_usage():
    """
    Basic usage example / 基本用法示例
    """
    print("=" * 60)
    print("GaussBio3D Example: DTI Feature Computation")
    print("GaussBio3D 示例: DTI特征计算")
    print("=" * 60)
    
    # NOTE: Replace these paths with your actual data files
    # 注意：将这些路径替换为您的实际数据文件
    pdb_path = "examples/target.pdb"
    sdf_path = "examples/drug.sdf"
    
    # Or use SMILES instead / 或使用SMILES代替
    # smiles = "CC(C)Cc1ccc(cc1)C(C)C(O)=O"  # Ibuprofen / 布洛芬
    
    print("\n1. Loading molecules... / 加载分子...")
    try:
        # Load protein / 加载蛋白质
        protein = Protein.from_pdb(pdb_path, chain_id="A")
        print(f"   Loaded protein with {len(protein.nodes)} atoms")
        print(f"   加载了包含 {len(protein.nodes)} 个原子的蛋白质")
        
        # Load ligand / 加载配体
        ligand = Ligand.from_sdf(sdf_path)
        print(f"   Loaded ligand with {len(ligand.nodes)} atoms")
        print(f"   加载了包含 {len(ligand.nodes)} 个原子的配体")
    except Exception as e:
        print(f"   Error loading files: {e}")
        print(f"   加载文件出错: {e}")
        print("   Please provide valid PDB and SDF files.")
        print("   请提供有效的PDB和SDF文件。")
        return
    
    print("\n2. Configuring mGLI parameters... / 配置mGLI参数...")
    config = MgliConfig(
        distance_bins=[0.0, 3.0, 6.0, 10.0, 20.0],  # Distance scales / 距离尺度
        use_rbf=False,                               # Use hard bins / 使用硬分箱
        signed=False,                                # Use absolute GLI / 使用绝对GLI
        stats=["sum", "mean", "max", "min"],        # Statistics / 统计量
        group_mode_A="residue_class",                # Protein grouping / 蛋白质分组
        group_mode_B="element",                      # Ligand grouping / 配体分组
    )
    print("   Configuration set / 配置完成")
    
    print("\n3. Computing global mGLI descriptor... / 计算全局mGLI描述符...")
    global_feat = global_mgli_descriptor(protein, ligand, config)
    print(f"   Global feature shape: {global_feat.shape}")
    print(f"   全局特征形状: {global_feat.shape}")
    print(f"   Feature dimension: {global_feat.shape[0]}")
    print(f"   特征维度: {global_feat.shape[0]}")
    
    print("\n4. Computing node-level features... / 计算节点级特征...")
    prot_node_feat = node_mgli_features(protein, ligand, config)
    lig_node_feat = node_mgli_features(ligand, protein, config)
    print(f"   Protein node features: {prot_node_feat.shape}")
    print(f"   蛋白质节点特征: {prot_node_feat.shape}")
    print(f"   Ligand node features: {lig_node_feat.shape}")
    print(f"   配体节点特征: {lig_node_feat.shape}")
    
    print("\n5. Computing pairwise mGLI matrix... / 计算成对mGLI矩阵...")
    pairwise_mat = pairwise_mgli_matrix(protein, ligand, signed=False, agg="mean")
    print(f"   Pairwise matrix shape: {pairwise_mat.shape}")
    print(f"   成对矩阵形状: {pairwise_mat.shape}")
    print(f"   Matrix statistics / 矩阵统计:")
    print(f"     Mean GLI: {np.mean(pairwise_mat):.6f}")
    print(f"     平均GLI: {np.mean(pairwise_mat):.6f}")
    print(f"     Max GLI: {np.max(pairwise_mat):.6f}")
    print(f"     最大GLI: {np.max(pairwise_mat):.6f}")
    
    print("\n6. Using convenience function... / 使用便捷函数...")
    print("   (Alternative way to compute all features at once)")
    print("   （一次性计算所有特征的替代方法）")
    
    # This computes all features in one call / 这会一次性计算所有特征
    # features = compute_dti_features(
    #     pdb_path=pdb_path,
    #     sdf_path=sdf_path,
    #     chain_id="A",
    #     config=config,
    # )
    
    print("\n" + "=" * 60)
    print("Example completed successfully! / 示例成功完成！")
    print("=" * 60)


def example_with_smiles():
    """
    Example using SMILES string / 使用SMILES字符串的示例
    """
    print("\n" + "=" * 60)
    print("GaussBio3D Example: Using SMILES")
    print("GaussBio3D 示例: 使用SMILES")
    print("=" * 60)
    
    # Ibuprofen SMILES / 布洛芬SMILES
    smiles = "CC(C)Cc1ccc(cc1)C(C)C(O)=O"
    
    print(f"\n1. Creating ligand from SMILES... / 从SMILES创建配体...")
    print(f"   SMILES: {smiles}")
    
    try:
        ligand = Ligand.from_smiles(smiles)
        print(f"   Ligand created with {len(ligand.nodes)} atoms")
        print(f"   创建了包含 {len(ligand.nodes)} 个原子的配体")
        print(f"   Number of bonds: {len(ligand.curves)}")
        print(f"   键的数量: {len(ligand.curves)}")
    except Exception as e:
        print(f"   Error: {e}")
        print(f"   错误: {e}")
        return
    
    print("\n2. Molecule information / 分子信息:")
    elements = [node.element for node in ligand.nodes]
    from collections import Counter
    element_counts = Counter(elements)
    for elem, count in sorted(element_counts.items()):
        print(f"   {elem}: {count}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    # Run the basic example / 运行基本示例
    # example_basic_usage()
    
    # Run the SMILES example / 运行SMILES示例
    example_with_smiles()
    
    print("\nNote: Uncomment example_basic_usage() if you have PDB/SDF files.")
    print("注意: 如果您有PDB/SDF文件，请取消注释example_basic_usage()。")
