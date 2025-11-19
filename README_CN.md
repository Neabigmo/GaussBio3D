# GaussBio3D: 多尺度高斯链接积分库

一个用于**小分子、蛋白质和核酸**的基于**多尺度高斯链接积分(mGLI)**的3D拓扑描述符Python库。

本库旨在为生物分子交互任务提供**统一的3D表示框架**，支持以下任务：

- 药物-靶点交互 (DTI)
- 蛋白质-蛋白质交互 (PPI)
- 药物-药物交互 (DDI)
- miRNA/核酸-靶点交互 (MTI)
- 蛋白质-DNA/RNA复合物等

---

## 1. 数学背景

### 1.1 高斯链接积分（连续形式）

给定两条光滑空间曲线 C₁ 和 C₂，**高斯链接积分**定义为：

```
GLI(C₁, C₂) = (1/4π) ∫∫ [(dr₁ × dr₂) · (r₁ - r₂)] / ||r₁ - r₂||³
              C₁ C₂
```

它度量两条曲线之间的**拓扑缠绕/缠结**关系。对于闭合曲线，它是一个整数（链接数），但对于开放曲线（如生物分子片段），它是一个实值的"链接强度"。

### 1.2 离散线段近似

我们用一组直线段来近似每条曲线：

- C₁ = {Lᵢ}, 其中 Lᵢ = [a₀, a₁]
- C₂ = {Mⱼ}, 其中 Mⱼ = [b₀, b₁]

则有：

```
GLI(C₁, C₂) ≈ Σᵢⱼ GLI(Lᵢ, Mⱼ)
```

对于线段 L=[a₀,a₁] 和 M=[b₀,b₁]，我们使用基于**球面几何的标准近似方法**。

---

## 2. 多尺度与分组mGLI特征

我们希望捕获分子A和B的各部分在**何种强度和何种距离尺度**下的拓扑链接特征。

### 2.1 节点对量

对于节点（原子/残基/碱基）i ∈ A, j ∈ B：

- 位置: xᵢ, xⱼ
- 距离: rᵢⱼ = ||xᵢ - xⱼ||
- 局部GLI: gᵢⱼ = 节点i和节点j相关联线段之间的聚合GLI

### 2.2 径向加权（多尺度）

我们定义径向基函数 φₖ(r)：

- **硬分箱**:
```
φₖ(r) = 𝟙[r ∈ [Rₖ, Rₖ₊₁)], k=1..K
```

- **RBF**:
```
φₖ(r) = exp(-(r-μₖ)²/(2σₖ²))
```

则多尺度聚合特征为：

```
hₖ = Σᵢⱼ φₖ(rᵢⱼ) · f(gᵢⱼ)
```

### 2.3 分组：元素/残基/碱基

我们进一步按离散类别对节点分组：

- 小分子: 元素/官能团
- 蛋白质: 残基类型或残基类别（疏水/芳香等）
- 核酸: 碱基类型(A/C/G/T/U)或主链vs碱基

---

## 3. 统一几何表示

我们将每个生物分子表示为：

- `Node` (节点): 原子/残基/碱基
- `Segment` (线段): 两个3D点之间的有向线段
- `Curve` (曲线): 由线段组成的折线
- `Structure` (结构): 节点+曲线的集合

---

## 4. 安装和依赖

需要：

- Python 3.9+
- `numpy`
- `rdkit` (用于SDF/MOL2/SMILES解析)
- `biopython` (用于PDB/mmCIF解析)

安装：

```bash
pip install numpy biopython rdkit-pypi
```

或从源码安装：

```bash
git clone https://github.com/yourusername/GaussBio3D
cd GaussBio3D
pip install -e .
```

---

## 5. 基本用法

### 5.1 计算蛋白质-配体全局mGLI描述符

```python
from gaussbio3d.molecules import Protein, Ligand
from gaussbio3d.config import MgliConfig
from gaussbio3d.features.descriptor import global_mgli_descriptor

# 加载蛋白质和配体
prot = Protein.from_pdb("examples/target.pdb", chain_id="A")
lig = Ligand.from_sdf("examples/drug.sdf")

# 配置mGLI参数
config = MgliConfig(
    distance_bins=[0.0, 3.0, 6.0, 10.0, 20.0],
    use_rbf=False,
    signed=False,
    group_mode_A="residue_class",
    group_mode_B="element",
)

# 计算全局描述符
feat = global_mgli_descriptor(prot, lig, config)
print("特征形状:", feat.shape)
```

### 5.2 DTI模型的节点级mGLI特征

```python
from gaussbio3d.features.node_features import node_mgli_features

# 计算节点级特征
node_feat_prot = node_mgli_features(prot, lig, config)
node_feat_lig  = node_mgli_features(lig, prot, config)
```

### 5.3 用于交叉注意力的成对mGLI矩阵

```python
from gaussbio3d.features.pairwise import pairwise_mgli_matrix

# 计算成对矩阵
M = pairwise_mgli_matrix(prot, lig, config)
# M.shape = (N_prot_nodes, N_lig_nodes)
```

---

## 6. 任务辅助工具

### DTI (药物-靶点交互)

```python
from gaussbio3d.tasks.dti import compute_dti_features

dti_feats = compute_dti_features(
    pdb_path="examples/target.pdb",
    sdf_path="examples/drug.sdf",
)
```

### PPI (蛋白质-蛋白质交互)

```python
from gaussbio3d.tasks.ppi import compute_ppi_features

ppi_feats = compute_ppi_features(
    pdb_path_A="protein_A.pdb",
    pdb_path_B="protein_B.pdb",
)
```

### MTI (分子-靶点交互)

```python
from gaussbio3d.tasks.mti import compute_mti_features

mti_feats = compute_mti_features(
    protein_pdb="protein.pdb",
    na_pdb="dna.pdb",
)
```

---

## 7. 项目结构

```
GaussBio3D/
├── gaussbio3d/
│   ├── __init__.py
│   ├── config.py              # 配置
│   ├── core/                  # 核心算法
│   │   ├── geometry.py        # 几何基元
│   │   └── gli.py             # GLI计算
│   ├── features/              # 特征提取
│   │   ├── descriptor.py      # 全局描述符
│   │   ├── node_features.py   # 节点级特征
│   │   └── pairwise.py        # 成对特征
│   ├── io/                    # 输入输出
│   │   ├── mol.py             # 分子文件I/O
│   │   └── pdb.py             # PDB文件I/O
│   ├── molecules/             # 分子表示
│   │   ├── ligand.py          # 小分子
│   │   ├── protein.py         # 蛋白质
│   │   └── nucleic_acid.py    # 核酸
│   └── tasks/                 # 特定任务辅助
│       ├── dti.py             # 药物-靶点交互
│       ├── ppi.py             # 蛋白质-蛋白质交互
│       └── mti.py             # 分子-靶点交互
├── examples/                  # 示例脚本
├── tests/                     # 单元测试
├── README.md
├── README_CN.md               # 中文说明
├── QUICKSTART.py              # 快速入门
├── setup.py
└── requirements.txt
```

---

## 8. 注意事项

* 本库是**研究原型**:
  * 效率尚未高度优化（GLI在最坏情况下是O(#segments²)）
  * 一些几何启发式方法被简化，应在生产使用中进一步优化

* 建议：
  * 根据您的任务调整距离分箱/RBF参数
  * 设计更细致的分组（如结合口袋残基vs非口袋残基）
  * 与因果/对抗训练流程集成以消除丰度偏差

---

## 9. 许可证

MIT License

---

## 10. 引用

如果您在研究中使用了GaussBio3D，请引用：

```bibtex
@software{gaussbio3d,
  title={GaussBio3D: Multiscale Gauss Linking Integral Library for Biomolecular 3D Topology},
  author={Your Name},
  year={2025},
  url={https://github.com/yourusername/GaussBio3D}
}
```

---

## 11. 联系方式

如有问题或建议，请通过以下方式联系：

- GitHub Issues: https://github.com/yourusername/GaussBio3D/issues
- Email: your.email@example.com
