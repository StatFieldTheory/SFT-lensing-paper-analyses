import sys, time, numpy as np
sys.path.insert(0, "..")
import _bootstrap
ds, bg_mod, core = _bootstrap.wire()
import corr_op_C_callable as co
tab = co._TABLE
lam = np.asarray(tab.cl_table.lambda_grid, float)
print("lambda_grid:", lam.size, "nodes", lam[:4], "...", lam[-4:])
n1 = np.array([0.,0.,1.])
g = 1.0/60.0*np.pi/180.0
n2 = np.array([np.sin(g),0.,np.cos(g)])
t0=time.time()
C = co.C_fn(n1, 1000.0, n2, 1500.0)
print("C(1000,1500) 00 =", C[0,0], " time", time.time()-t0)
C2 = co.C_fn(n2, 1500.0, n1, 1000.0)
print("swap-check 00 =", C2[0,0], "rel", abs(C[0,0]-C2[0,0])/abs(C[0,0]))
C3 = co.C_fn(n1, 1500.0, n2, 1000.0)
print("t-swap same rays 00 =", C3[0,0], "rel vs C", abs(C[0,0]-C3[0,0])/abs(C[0,0]))
# batch timing
NN=64
t1 = np.repeat(np.linspace(406,2313,NN), NN)
t2 = np.tile(np.linspace(406,2313,NN), NN)
t0=time.time()
B = co.C_fn_batch(n1, t1, n2, t2)
print("batch", B.shape, "time", time.time()-t0)
