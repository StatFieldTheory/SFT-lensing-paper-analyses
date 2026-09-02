"""Read the paper's converged permutation-closed FK fold (order 2, a=b=0)."""
import numpy as np

PATH = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/"
        "sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz")


def fk_table():
    d = np.load(PATH, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    xs, ys, vs = d["x"][m], d["y"][m], d["value"][m]
    out = {}
    for xi, yi, v in zip(xs, ys, vs):
        c = float(np.clip(np.dot(np.asarray(xi, float), np.asarray(yi, float)), -1, 1))
        gam = float(np.degrees(np.arccos(c)) * 60.0)
        out[round(gam, 4)] = float(v)
    return dict(sorted(out.items()))


if __name__ == "__main__":
    for g, v in fk_table().items():
        print(f"{g:10.4f}'  {v: .6e}")
