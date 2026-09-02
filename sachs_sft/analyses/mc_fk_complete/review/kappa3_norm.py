"""Exact third moment of the LINEAR observable -- no F vertex anywhere.

kappa1_A = sum_j dlam Wd[j] f_A[j],   f = z + 0.5 Q:(zz - C)

At O(Q) the connected third moment is EXACTLY

    <k1_A k1_B k1_C>_c = sum_j dlam Wd_j [ Q[j]_Amn Y[j]_mB Y[j]_nC + (A<->B) + (A<->C) ]
    Y[j]_mb = sum_p dlam Wd_p R[j,p]_mb ,     R[j,p] = rho^|j-p| A[min(j,p)]

and the local (white-noise) prediction it must reproduce is

    sum_j dlam Wd_j^3 zeta[j].

If Wd varies slowly over a correlation length, Y[j] -> Wd_j B[j] with
B[j] = sum_p dlam R[j,p], so the requirement collapses to the node-wise
condition  cum3_from_Q(Q[j], B[j]) = zeta[j],  i.e. Q = solve_Q(B, zeta).

Everything below is O(N) in memory: Y and B are built by two-sided recursions
rather than by forming the (N,N,6,6) correlation.
"""
import numpy as np


def two_sided(A, rho, dlam, weight):
    """S[j] = sum_p dlam weight[p] rho^|j-p| A[min(j,p)]   -- exact, O(N)."""
    N = A.shape[0]
    left = np.zeros_like(A)                    # sum_{p<j}
    for j in range(1, N):
        left[j] = rho * (left[j - 1] + dlam * weight[j - 1] * A[j - 1])
    c2 = np.empty(N)                           # sum_{p>=j} dlam w_p rho^{p-j}
    c2[N - 1] = dlam * weight[N - 1]
    for j in range(N - 2, -1, -1):
        c2[j] = dlam * weight[j] + rho * c2[j + 1]
    return left + c2[:, None, None] * A


def kappa1_third_moment(grid, A, Q, rho):
    """Exact O(Q) third cumulant of the linear observable, (6,6,6)."""
    Y = two_sided(A, rho, grid.dlam, grid.Wd)
    t1 = np.einsum("j,jamn,jmb,jnc->abc", grid.dlam * grid.Wd, Q, Y, Y,
                   optimize=True)
    return t1 + t1.transpose(1, 0, 2) + t1.transpose(2, 1, 0)


def local_prediction(grid, zeta):
    return np.einsum("j,jabc->abc", grid.dlam * grid.Wd ** 3, zeta, optimize=True)


def smeared_B(grid, A, rho):
    """Same as fk_expect_exact.smeared_covariance, but O(N) instead of O(N^2)."""
    return two_sided(A, rho, grid.dlam, np.ones(grid.lam.size))
