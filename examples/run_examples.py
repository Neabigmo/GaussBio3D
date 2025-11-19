"""
Run feature extraction on provided example files and save results.
在提供的示例文件上运行特征提取并保存结果。

This script computes DTI features for:
- Protein (mmCIF) + Ligand (SDF)
- Protein (PDB)   + Ligand (SDF)

Output files are saved under `examples/result`.
"""

import os
import json
import numpy as np

from gaussbio3d.tasks.dti import compute_dti_features
from gaussbio3d.config import MgliConfig


def main():
    # Absolute paths provided by user
    sdf_path = r"F:\1-2025project\GAC\GaussBio3D\examples\0beaa2fb7ef0.sdf"
    cif_path = r"F:\1-2025project\GAC\GaussBio3D\examples\1AGW-assembly1.cif"
    pdb_path = r"F:\1-2025project\GAC\GaussBio3D\examples\AF-Q5TBG4-F1-model_v6.pdb"
    out_dir = r"F:\1-2025project\GAC\GaussBio3D\examples\result"

    os.makedirs(out_dir, exist_ok=True)

    config = MgliConfig()

    # 1) CIF + SDF
    dti_cif = compute_dti_features(
        pdb_path=cif_path,
        sdf_path=sdf_path,
        chain_id=None,
        config=config,
    )
    np.savez_compressed(
        os.path.join(out_dir, "dti_cif_sdf.npz"),
        global_feat=dti_cif["global_feat"],
        prot_node_feat=dti_cif["prot_node_feat"],
        lig_node_feat=dti_cif["lig_node_feat"],
        pairwise_mgli=dti_cif["pairwise_mgli"],
    )

    # 2) PDB + SDF
    dti_pdb = compute_dti_features(
        pdb_path=pdb_path,
        sdf_path=sdf_path,
        chain_id="A",  # the provided PDB uses chain A
        config=config,
    )
    np.savez_compressed(
        os.path.join(out_dir, "dti_pdb_sdf.npz"),
        global_feat=dti_pdb["global_feat"],
        prot_node_feat=dti_pdb["prot_node_feat"],
        lig_node_feat=dti_pdb["lig_node_feat"],
        pairwise_mgli=dti_pdb["pairwise_mgli"],
    )

    # Write a small JSON summary for quick inspection
    summary = {
        "dti_cif_sdf": {
            "global_feat_shape": list(dti_cif["global_feat"].shape),
            "prot_node_feat_shape": list(dti_cif["prot_node_feat"].shape),
            "lig_node_feat_shape": list(dti_cif["lig_node_feat"].shape),
            "pairwise_mgli_shape": list(dti_cif["pairwise_mgli"].shape),
        },
        "dti_pdb_sdf": {
            "global_feat_shape": list(dti_pdb["global_feat"].shape),
            "prot_node_feat_shape": list(dti_pdb["prot_node_feat"].shape),
            "lig_node_feat_shape": list(dti_pdb["lig_node_feat"].shape),
            "pairwise_mgli_shape": list(dti_pdb["pairwise_mgli"].shape),
        },
    }
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)

    print("Results saved to:", out_dir)


if __name__ == "__main__":
    main()