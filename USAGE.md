# GaussBio3D Library Usage Guide / 使用指南

This guide shows how to call GaussBio3D via Python code: loading structures, configuring mGLI, computing descriptors, and using presets/pipeline.

本文档介绍如何通过 Python 代码调用 GaussBio3D：加载结构、配置 mGLI、计算特征，以及使用预设与流水线。

## Install / 安装

```bash
pip install gaussbio3d
# 或 Conda 安装 RDKit 后再安装包
conda install -c conda-forge rdkit
pip install gaussbio3d
```

## Load Structures / 加载结构

```python
from gaussbio3d.molecules import Protein, Ligand

prot = Protein.from_pdb("protein.pdb", chain_id="A")
lig  = Ligand.from_sdf("ligand.sdf")
```

## Configure mGLI / 配置 mGLI

```python
from gaussbio3d.config import MgliConfig

cfg = MgliConfig(
    distance_bins=[3.0, 5.0, 7.0, 10.0, 15.0],
    use_rbf=False,
    signed=False,
    group_mode_A="element",
    group_mode_B="element",
    max_distance=15.0,
    n_jobs=4,
)
```

## Global Descriptor / 全局描述符

```python
from gaussbio3d.features.descriptor import global_mgli_descriptor

feat = global_mgli_descriptor(prot, lig, cfg)
print(feat.shape)  # 扁平向量维度
```

## Pairwise Matrix / 成对矩阵

```python
from gaussbio3d.features.pairwise import pairwise_mgli_matrix

M = pairwise_mgli_matrix(prot, lig, cfg)
print(M.shape)  # (N_prot_nodes, N_lig_nodes)
```

## Node Features / 节点级特征

```python
from gaussbio3d.features.node_features import node_mgli_features

node_feat_prot = node_mgli_features(prot, lig, cfg)
node_feat_lig  = node_mgli_features(lig, prot, cfg)
```

## Presets / 预设入口

```python
from gaussbio3d import (
    flexibility_mgli_pipeline,
    compute_bfactor_mgli,
    default_dti_mgli_pipeline,
    compute_pl_complex_mgli,
)

# 一键计算蛋白–配体互作特征
vec_pl = compute_pl_complex_mgli("protein.pdb", "ligand.sdf")

# 一键计算蛋白柔性特征
vec_flex = compute_bfactor_mgli("protein.cif")
```

## Pipeline / 流水线

```python
from gaussbio3d.core.pipeline import MGLIPipeline, PCAProjector

pipe = MGLIPipeline(config=cfg, projector=PCAProjector(n_components=256))
X = pipe.fit_transform([(prot, lig)])
print(X.shape)
```

