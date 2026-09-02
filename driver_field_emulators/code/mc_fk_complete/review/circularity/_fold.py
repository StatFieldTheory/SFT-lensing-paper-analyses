"""Read an sft-wick FK fold npz at its OWN gamma nodes (no interpolation)."""
import math
import numpy as np


def fold_nodes(path, t_final=None, a=0, b=0, order=2):
    d = np.load(path, allow_pickle=True)
    m = (d["a"] == a) & (d["b"] == b) & (d["order"] == order)
    if t_final is not None:
        m &= np.isclose(d["t_final"], t_final)
    g = []
    for x, y in zip(d["x"][m], d["y"][m]):
        x = np.asarray(x, float); y = np.asarray(y, float)
        g.append(math.degrees(math.acos(float(np.clip(
            np.dot(x / np.linalg.norm(x), y / np.linalg.norm(y)), -1, 1)))) * 60)
    g = np.asarray(g); v = np.asarray(d["value"], float)[m]
    o = np.argsort(g)
    return g[o], v[o]
