"""Feynman-diagram figures of the manuscript: the Order-2 two-point diagrams
(FF, FK, KK) and the Order-1 three-point diagrams.

Provenance. The deployed ``figures/FeynDiag_2pt_order2.pdf`` and
``figures/FeynDiag_3pt_order1.pdf`` are byte-identical to
``sft-wick/examples/FeynDiag_*.pdf`` (md5 ``c26347c4...`` and ``66a713ae...``),
which cells 29 and 38 of ``sft-wick/examples/nonlocal_vertex_2pt.ipynb`` wrote
on 2026-06-03 (sft-wick commit ``dd03748``). This script is that notebook's
generator path extracted verbatim (cells 2, 3, 5, 7, 24, 29, 33, 38), with the
diagnostic prints dropped and the output redirected. The rendering code of
sft-wick has no commits after ``dd03748`` (checked 2026-09-02).

Output (default): ``figures/outputs/FeynDiag_2pt_order2.pdf`` and
``figures/outputs/FeynDiag_3pt_order1.pdf`` next to this script; the paper
copies under the repository ``figures/`` are deployed by the reproduction
driver. Pass ``--out-dir`` to redirect.

Run (sft-wick env, sft_wick importable):
    /opt/homebrew/Caskroom/miniconda/base/envs/sft-wick/bin/python \
        SFT-lensing-paper-analyses/figures/make_feyndiag.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from sft_wick import Action, Field, Vertex, compute_moment, reset_uid_counter  # noqa: E402

HERE = Path(__file__).resolve().parent


def mathcal_x_external_label(node_id, attrs):
    """Notebook-specific label override: physical external fields as calligraphic X."""
    field_type = attrs.get("field_type")
    field_type = field_type.value if hasattr(field_type, "value") else field_type
    if field_type != "physical":
        return None
    comp = attrs.get("component")
    if comp is None:
        return r"$\mathcal{X}$"
    return rf"$\mathcal{{X}}_{{{comp}}}$"


def build_action():
    reset_uid_counter()  # reproducible FieldOperator UIDs, as in the notebook
    X = Field("X", "physical", n_components=3)
    psi = Field("psi", "response", n_components=3)
    # Vertex 1: F_{abc} psi_a X_b X_c (local)
    v1 = Vertex(fields=[psi, X, X], coupling="F")
    # Vertex 2: K_{abc}(x1,x2,x3) psi_a(x1) psi_b(x2) psi_c(x3) (non-local)
    v2 = Vertex(fields=[psi, psi, psi], coupling="K", local=False)
    return X, Action(vertices=[v1, v2])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=HERE / "outputs")
    a = ap.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)

    X, action = build_action()

    # Two-point correlator at order 2 (cells 5, 7, 29).
    obs_2pt = [X("a", "x"), X("b", "y")]
    result = compute_moment(obs_2pt, action, order=2, response_phase=True)
    fig1 = result.draw_diagrams(order=2, figsize=(4, 3.0),
                                external_label_fn=mathcal_x_external_label, ncols=3)
    out1 = a.out_dir / "FeynDiag_2pt_order2.pdf"
    fig1.savefig(out1, bbox_inches="tight")
    print(f"wrote {out1}")

    # Three-point correlator at order 1 (cells 33, 38).
    obs_3pt = [X("a", "x"), X("b", "y"), X("c", "z")]
    result_3pt = compute_moment(obs_3pt, action, order=1, response_phase=True)
    fig2 = result_3pt.draw_diagrams(order=1, figsize=(4, 3), ncols=3,
                                    external_label_fn=mathcal_x_external_label)
    out2 = a.out_dir / "FeynDiag_3pt_order1.pdf"
    fig2.savefig(out2, bbox_inches="tight")
    print(f"wrote {out2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
