"""
Molecule file I/O using RDKit
使用RDKit的分子文件输入/输出

This module provides functions to load small molecules from various file formats
(SDF, MOL2, SMILES) and extract coordinates and connectivity information.

本模块提供从各种文件格式(SDF、MOL2、SMILES)加载小分子并提取坐标和连接信息的函数。
"""

from __future__ import annotations

from typing import Tuple, List
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:  # pragma: no cover
    Chem = None
    AllChem = None


def load_mol_from_sdf(path: str) -> "Chem.Mol":
    """
    Load an RDKit mol object from an SDF file.
    从SDF文件加载RDKit分子对象。
    
    Parameters / 参数
    ----------
    path : str
        Path to SDF file / SDF文件路径
        
    Returns / 返回
    -------
    Chem.Mol
        RDKit molecule object / RDKit分子对象
        
    Raises / 引发
    ------
    ImportError
        If RDKit is not installed / 如果未安装RDKit
    ValueError
        If no valid molecule found in file / 如果文件中未找到有效分子
    """
    if Chem is None:
        raise ImportError("RDKit is required for SDF parsing (pip install rdkit-pypi).")
    suppl = Chem.SDMolSupplier(path, removeHs=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        raise ValueError(f"No valid molecule found in SDF: {path}")
    return mols[0]


def simple_parse_sdf_v2000(path: str) -> Tuple[np.ndarray, List[str], List[Tuple[int, int]]]:
    """
    Minimal SDF V2000 parser to extract coordinates, elements, and bonds.
    最小化的SDF V2000解析器，用于提取坐标、元素和键。

    This function is a fallback when RDKit is unavailable.
    当RDKit不可用时，本函数作为后备方案。

    Parameters / 参数
    ----------
    path : str
        Path to SDF file / SDF文件路径

    Returns / 返回
    -------
    coords : np.ndarray
        Atomic coordinates array of shape (N, 3)
    elements : List[str]
        Element symbol per atom
    bonds : List[Tuple[int, int]]
        Zero-based atom index pairs for bonds
    """
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        lines = [ln.rstrip("\n") for ln in f]

    # Find counts line (typically 4th line) containing V2000
    counts_idx = None
    for i in range(min(10, len(lines))):
        if "V2000" in lines[i]:
            counts_idx = i
            break
    if counts_idx is None:
        # Fallback: search globally
        for i, ln in enumerate(lines):
            if "V2000" in ln:
                counts_idx = i
                break
    if counts_idx is None:
        raise ValueError(f"SDF V2000 counts line not found: {path}")

    parts = lines[counts_idx].split()
    if len(parts) < 2:
        raise ValueError(f"Malformed counts line in SDF: {path}")
    try:
        num_atoms = int(parts[0])
        num_bonds = int(parts[1])
    except Exception as e:
        raise ValueError(f"Failed to parse atom/bond counts in SDF: {path}") from e

    atom_start = counts_idx + 1
    bond_start = atom_start + num_atoms

    coords = np.zeros((num_atoms, 3), dtype=float)
    elements: List[str] = []
    for i in range(num_atoms):
        ln = lines[atom_start + i]
        # x y z symbol ... fixed-width or space-separated
        try:
            parts = ln.split()
            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
            elem = parts[3]
        except Exception as e:
            raise ValueError(f"Malformed atom line at {atom_start+i} in SDF: {path}") from e
        coords[i] = [x, y, z]
        elements.append(elem)

    bonds: List[Tuple[int, int]] = []
    for i in range(num_bonds):
        ln = lines[bond_start + i]
        parts = ln.split()
        if len(parts) < 2:
            # Some SDFs may have trailing blocks; stop when hit non-bond lines
            break
        try:
            a = int(parts[0]) - 1
            b = int(parts[1]) - 1
        except Exception:
            break
        if 0 <= a < num_atoms and 0 <= b < num_atoms:
            bonds.append((a, b))

    return coords, elements, bonds


def load_mol_from_mol2(path: str) -> "Chem.Mol":
    """
    Load an RDKit mol object from a MOL2 file.
    从MOL2文件加载RDKit分子对象。
    
    Parameters / 参数
    ----------
    path : str
        Path to MOL2 file / MOL2文件路径
        
    Returns / 返回
    -------
    Chem.Mol
        RDKit molecule object / RDKit分子对象
        
    Raises / 引发
    ------
    ImportError
        If RDKit is not installed / 如果未安装RDKit
    ValueError
        If failed to read MOL2 file / 如果读取MOL2文件失败
    """
    if Chem is None:
        raise ImportError("RDKit is required for MOL2 parsing.")
    mol = Chem.MolFromMol2File(path, removeHs=False)
    if mol is None:
        raise ValueError(f"Failed to read MOL2: {path}")
    return mol


def load_mol_from_smiles(smiles: str) -> "Chem.Mol":
    """
    Load an RDKit mol object from SMILES and generate 3D conformer.
    从SMILES加载RDKit分子对象并生成3D构象。
    
    Parameters / 参数
    ----------
    smiles : str
        SMILES string / SMILES字符串
        
    Returns / 返回
    -------
    Chem.Mol
        RDKit molecule object with 3D coordinates / 带3D坐标的RDKit分子对象
        
    Raises / 引发
    ------
    ImportError
        If RDKit is not installed / 如果未安装RDKit
    ValueError
        If failed to parse SMILES / 如果解析SMILES失败
    """
    if Chem is None or AllChem is None:
        raise ImportError("RDKit is required for SMILES/3D generation.")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def mol_to_coordinates_and_elements(mol: "Chem.Mol") -> Tuple[np.ndarray, List[str]]:
    """
    Extract 3D coordinates and element symbols from an RDKit molecule.
    从RDKit分子中提取3D坐标和元素符号。

    Parameters / 参数
    ----------
    mol : Chem.Mol
        RDKit molecule object / RDKit分子对象

    Returns / 返回
    -------
    coords : np.ndarray
        Atomic coordinates, shape (N_atoms, 3) / 原子坐标，形状为(N_atoms, 3)
    elements : List[str]
        Element symbol per atom / 每个原子的元素符号
    """
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    coords = np.zeros((n, 3), dtype=float)
    elements: List[str] = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        coords[i] = [pos.x, pos.y, pos.z]
        elements.append(atom.GetSymbol())
    return coords, elements


def mol_to_bond_pairs(mol: "Chem.Mol") -> List[Tuple[int, int]]:
    """
    Extract bond pairs (i, j) from RDKit molecule.
    从RDKit分子中提取键对(i, j)。
    
    Parameters / 参数
    ----------
    mol : Chem.Mol
        RDKit molecule object / RDKit分子对象
        
    Returns / 返回
    -------
    bonds : List[Tuple[int, int]]
        List of atom index pairs connected by bonds
        由键连接的原子索引对列表
    """
    bonds = []
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        bonds.append((i, j))
    return bonds
