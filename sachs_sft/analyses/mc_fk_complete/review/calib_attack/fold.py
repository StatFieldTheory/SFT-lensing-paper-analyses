import numpy as np
PATH = ("/Users/zzhang/Documents/MyDrafts/STF_lensing/driver_field_emulators/"
        "sftwick_outputs/2PCF/C_corr_op_K_limber_FK_cut15360_permfix/xi_C_corr_op_K_limber_FK_cut15360_permfix.npz")

def fk_table():
    d = np.load(PATH, allow_pickle=True)
    m = (d["a"] == 0) & (d["b"] == 0) & (d["order"] == 2)
    out = {}
    for xi, yi, v in zip(d["x"][m], d["y"][m], d["value"][m]):
        c = float(np.clip(np.dot(np.asarray(xi, float), np.asarray(yi, float)), -1, 1))
        out[float(np.degrees(np.arccos(c)) * 60.0)] = float(v)
    return dict(sorted(out.items()))

def fold_at(gam):
    t = fk_table()
    gs = np.array(list(t.keys())); vs = np.array(list(t.values()))
    # exact node if available
    i = int(np.argmin(np.abs(gs - gam)))
    if abs(gs[i] - gam) < 1e-3:
        return vs[i], True
    return float(np.exp(np.interp(np.log(gam), np.log(gs), np.log(np.abs(vs))))), False

if __name__ == "__main__":
    for g, v in fk_table().items():
        print(f"{g:10.4f}'  {v: .6e}")
