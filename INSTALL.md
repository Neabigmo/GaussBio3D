# Installation Guide / 安装指南

## Method 1: Install from Source (Recommended) / 方法1：从源码安装（推荐）

### Step 1: Clone or Download / 步骤1：克隆或下载

```bash
# If using git / 如果使用git
git clone https://github.com/yourusername/GaussBio3D.git
cd GaussBio3D

# Or download and extract the zip file / 或下载并解压zip文件
```

### Step 2: Install Dependencies / 步骤2：安装依赖

```bash
# Install required packages / 安装必需的包
pip install numpy biopython rdkit-pypi

# Or install from requirements.txt / 或从requirements.txt安装
pip install -r requirements.txt
```

### Step 3: Install GaussBio3D / 步骤3：安装GaussBio3D

```bash
# Install in development mode (editable) / 以开发模式安装（可编辑）
pip install -e .

# Or install normally / 或正常安装
pip install .
```

---

## Method 2: Direct Installation / 方法2：直接安装

```bash
pip install numpy biopython rdkit-pypi
pip install git+https://github.com/yourusername/GaussBio3D.git
```

---

## Verify Installation / 验证安装

### Quick Test / 快速测试

```python
# Open Python interpreter / 打开Python解释器
python

# Try importing / 尝试导入
>>> import gaussbio3d
>>> print(gaussbio3d.__version__)
0.1.0

>>> from gaussbio3d.molecules import Protein, Ligand
>>> from gaussbio3d.config import MgliConfig
>>> print("Installation successful! / 安装成功！")
```

### Run Example / 运行示例

```bash
# Run the SMILES example / 运行SMILES示例
python examples/example_dti.py
```

---

## Troubleshooting / 故障排除

### Issue 1: RDKit Installation Fails / 问题1：RDKit安装失败

**Solution / 解决方案:**

```bash
# Try conda installation / 尝试conda安装
conda install -c conda-forge rdkit

# Or try different version / 或尝试不同版本
pip install rdkit-pypi==2022.9.5
```

### Issue 2: Biopython Import Error / 问题2：Biopython导入错误

**Solution / 解决方案:**

```bash
pip install --upgrade biopython
```

### Issue 3: NumPy Version Conflict / 问题3：NumPy版本冲突

**Solution / 解决方案:**

```bash
pip install --upgrade numpy
```

---

## System Requirements / 系统要求

- **Python**: 3.9 or higher / 3.9或更高版本
- **Operating System / 操作系统**: 
  - Windows 10/11
  - macOS 10.14+
  - Linux (Ubuntu 18.04+, CentOS 7+)
- **Memory / 内存**: Minimum 4GB RAM / 最少4GB内存
- **Disk Space / 磁盘空间**: ~500MB (including dependencies / 包括依赖项)

---

## Optional Dependencies / 可选依赖

For development and testing / 用于开发和测试:

```bash
pip install pytest pytest-cov black flake8
```

For Jupyter notebooks / 用于Jupyter notebooks:

```bash
pip install jupyter ipython

### GPU & Topology (Optional) / GPU与拓扑（可选）

To enable GPU acceleration (PyTorch) and PH topology features, install optional packages / 如需启用GPU加速（PyTorch）与PH拓扑特性，请安装可选包：

```bash
# Install optional dependencies individually / 单独安装
pip install numba ripser

# Install PyTorch (choose CUDA or CPU per your environment) / 安装PyTorch（按环境选择CUDA或CPU）
# CUDA 12.1 example / CUDA 12.1示例：
pip install torch --index-url https://download.pytorch.org/whl/cu121

# Or CPU-only / 或仅CPU：
pip install torch --index-url https://download.pytorch.org/whl/cpu

# Alternatively, install from optional requirements file / 或通过可选依赖文件安装
pip install -r requirements-optional.txt
```

After installation, set in config to use these features / 安装后在配置中启用：

```python
from gaussbio3d.config import MgliConfig
cfg = MgliConfig(max_distance=8.0, n_jobs=4, use_gpu=True)
```

Optional packages summary / 可选包摘要：
- `numba`：JIT优化部分CPU路径（可选）
- `torch`：启用GPU后端（可选）
- `ripser`：持久同调（PH）特性（可选）
```

---

## Uninstallation / 卸载

```bash
pip uninstall gaussbio3d
```

---

## Getting Help / 获取帮助

If you encounter any issues / 如果遇到任何问题:

1. Check the documentation / 查看文档:
   - README.md (English)
   - README_CN.md (中文)
   - QUICKSTART.py

2. Search existing issues / 搜索现有问题:
   https://github.com/yourusername/GaussBio3D/issues

3. Create a new issue / 创建新问题:
   https://github.com/yourusername/GaussBio3D/issues/new

4. Contact / 联系:
   your.email@example.com
