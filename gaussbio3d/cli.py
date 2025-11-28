from __future__ import annotations

import argparse
import numpy as np

from .presets import (
    default_dti_mgli_pipeline,
    flexibility_mgli_pipeline,
    compute_pl_complex_mgli,
    compute_bfactor_mgli,
)
from .molecules.protein import Protein
from .molecules.ligand import Ligand


def _cmd_compute(args: argparse.Namespace) -> int:
    mode = args.mode.strip().lower()
    out = args.out
    if mode == "pl":
        if not args.protein or not args.ligand:
            raise SystemExit("--protein and --ligand are required for mode=pl")
        vec = compute_pl_complex_mgli(args.protein, args.ligand, chain_id=args.chain)
        print(f"descriptor shape: {vec.shape}")
        if out:
            np.save(out, vec)
            print(f"saved to: {out}")
        return 0
    elif mode in {"protein-flex", "flex"}:
        if not args.protein:
            raise SystemExit("--protein is required for mode=protein-flex")
        vec = compute_bfactor_mgli(args.protein, chain_id=args.chain)
        print(f"descriptor shape: {vec.shape}")
        if out:
            np.save(out, vec)
            print(f"saved to: {out}")
        return 0
    else:
        raise SystemExit(f"unknown mode: {mode}")


def main() -> int:
    parser = argparse.ArgumentParser(prog="gaussbio3d", description="GaussBio3D CLI")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_compute = sub.add_parser("compute", help="compute mGLI descriptors")
    p_compute.add_argument("--mode", required=True, help="pl | protein-flex")
    p_compute.add_argument("--protein", help="protein PDB/mmCIF path")
    p_compute.add_argument("--ligand", help="ligand SDF path")
    p_compute.add_argument("--chain", help="protein chain id")
    p_compute.add_argument("--out", help="output .npy path")
    p_compute.set_defaults(func=_cmd_compute)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())

