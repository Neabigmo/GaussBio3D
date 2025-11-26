# 文档结构
## 概览
- 介绍 GaussBio3D 的定位：基于多尺度高斯链接积分（mGLI）的三维拓扑特征库，统一表征蛋白、配体、核酸的空间交互，用于 DTI/PPI/MTI 等任务。
- 指出三层输出：全局描述符、节点级特征、节点成对矩阵。

## 模式与功能
- 径向模式：硬分箱与 RBF 两种（`gaussbio3d/config.py:12`、`gaussbio3d/features/descriptor.py:129`）。
- 分组模式：元素分组与生物学分组（残基类别/碱基类型），参见 `gaussbio3d/molecules/protein.py:21`、`gaussbio3d/molecules/nucleic_acid.py:41`。
- 符号模式：`signed` 控制 GLI 的手性（`gaussbio3d/core/gli_segment.py:212`）。
- 统计聚合：`sum/mean/max/min/median`。

## 架构与实现
- 几何与结构：`Structure`/`Segment`/`Curve`（`gaussbio3d/core/geometry.py`）。
- 线段级 GLI：Numba 加速与 Numpy 回退（`gaussbio3d/core/gli_segment.py:218`、`gaussbio3d/core/gli_segment.py:232`）。
- 成对节点 GLI：距离剪枝、并行/GPU 可选（`gaussbio3d/core/pairwise_gli.py:26`）。
- 特征层：全局（`gaussbio3d/features/descriptor.py:129`）、节点级（`gaussbio3d/features/node_features.py:20`）、成对矩阵（`gaussbio3d/features/pairwise.py:19`）。
- 任务封装：DTI（`gaussbio3d/tasks/dti.py:23`）等入口。

## 关键 API
- `global_mgli_descriptor`、`node_mgli_features`、`pairwise_mgli_matrix` 的用途、输入/输出维度与常用参数；配以简短代码片段。
- DTI 快速特征：`compute_dti_features`（加载 PDB/SDF/SMILES 并生成四类输出）。

## 配置项
- 列出 `MgliConfig` 字段、默认值与效果（`gaussbio3d/config.py:12-78`），包含序列化接口（`gaussbio3d/config.py:79-86`）。

## 使用示例
- 构建分子结构：`Protein.from_pdb`、`Ligand.from_sdf/from_smiles`、`NucleicAcid.from_pdb`。
- 计算三层特征与保存示例。

## 性能与依赖
- 依赖：`numpy`、`biopython`、`rdkit-pypi`、可选 `numba`、`torch`、`ripser`。
- 性能开关：`max_distance` 距离剪枝、`n_jobs` 并行、`use_gpu`。

## 扩展与示例
- `examples/run_examples.py` 的说明与结果格式。
- 拓展：PH 等拓扑增强模块的可选性。

# 交付方式
- 将内容写入中文文档 `README_CN.md` 的“深入说明”段，或新增 `docs/overview.md`（中文），并在 `README.md` 添加指向链接。
- 在 `examples/` 中补充一段最小可运行示例（若现有示例足够，仅在文档引用）。

# 验证
- 本地运行 `examples/run_examples.py` 验证生成特征维度与落盘结果。
- 逐文件检查代码引用的行号与函数签名一致性。

# 后续迭代
- 根据使用反馈补充 FAQ 与最佳实践（分箱 vs RBF、何时使用 `signed`、如何选择 `group_mode`）。

# 需要确认
- 文档落点：更新 `README_CN.md` 还是新增 `docs/overview.md`。
- 是否需要英文版同步（`README.md` 增补“Deep Dive”）。