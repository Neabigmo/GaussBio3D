"""
Quick synthetic check of GLI computations without RDKit/Biopython.
无需RDKit/Biopython的GLI快速合成校验。
"""

import numpy as np
from gaussbio3d.core.geometry import Node, Segment, Curve, Structure
from gaussbio3d.core.pairwise_gli import compute_pairwise_node_gli


def make_simple_structure(name: str) -> Structure:
    s = Structure(metadata={"name": name})
    # Two nodes
    s.add_node(Node(id=0, coord=np.array([0.0, 0.0, 0.0]), element="C", group="C"))
    s.add_node(Node(id=1, coord=np.array([1.0, 0.0, 0.0]), element="C", group="C"))
    # One segment between them as two half-bonds
    midpoint = 0.5 * (s.nodes[0].coord + s.nodes[1].coord)
    seg1 = Segment(start=s.nodes[0].coord, end=midpoint, start_node_id=0, end_node_id=None, start_type="C", end_type="C")
    seg2 = Segment(start=s.nodes[1].coord, end=midpoint, start_node_id=1, end_node_id=None, start_type="C", end_type="C")
    s.add_curve(Curve(segments=[seg1, seg2], curve_type="bond"))
    return s


def main():
    A = make_simple_structure("A")
    B = make_simple_structure("B")
    gij, rij = compute_pairwise_node_gli(A, B, signed=False, agg="mean")
    print("gij shape:", gij.shape)
    print("rij shape:", rij.shape)
    print("gij:\n", gij)
    print("rij:\n", rij)


if __name__ == "__main__":
    main()