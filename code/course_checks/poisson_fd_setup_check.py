"""Environment sanity check for the BVP stack (NOT a coding task).

Distributed as-is: no TODOs to fill.  It verifies that the sparse
assembly + spsolve path used across Project 2 works on your machine
BEFORE you start the Lecture 3 coding lab, and prints the algebraic
residual that the Environment Setup addendum asks you to record.

Run:
    python poisson_fd_setup_check.py

Expected output (n = 15):
    algebraic residual ||b - A u_h||_inf = ~1e-14   (near machine precision)
    centre u_h(0.5, 0.5) = ~1.003                    (PROBLEM = P1 sanity value)

This solves -Delta u = 2 pi^2 sin(pi x) sin(pi y) with zero Dirichlet
data, whose exact solution is u* = sin(pi x) sin(pi y) (u*(0.5, 0.5) = 1).
The same five-point matrix, RHS convention, and PROBLEM = P1 case appear
in poisson_fd_skeleton.py --- running this script proves that code path
works before you add your assembly to the scaffold.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def assemble_poisson(n):
    """Five-point -Delta on (0, 1)^2, zero Dirichlet, b = f at interior nodes.

    Mirrors poisson_fd_skeleton.py: k = i*n + j with meshgrid indexing='ij'
    (axis 0 is x), f = 2 pi^2 sin(pi x) sin(pi y) so u* = sin(pi x) sin(pi y)
    vanishes on the whole boundary and b = f needs no boundary-value shift.
    """
    h = 1.0 / (n + 1)
    grid = np.arange(1, n + 1) * h
    X, Y = np.meshgrid(grid, grid, indexing="ij")
    rows, cols, data = [], [], []
    b = (2.0 * np.pi**2 * np.sin(np.pi * X) * np.sin(np.pi * Y)).ravel()
    k = lambda i, j: i * n + j
    for i in range(n):
        for j in range(n):
            rows.append(k(i, j))
            cols.append(k(i, j))
            data.append(4.0 / h**2)
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + di, j + dj
                if 0 <= ii < n and 0 <= jj < n:
                    rows.append(k(i, j))
                    cols.append(k(ii, jj))
                    data.append(-1.0 / h**2)
    A = coo_matrix((data, (rows, cols)), shape=(n * n, n * n)).tocsr()
    return A, b, h


if __name__ == "__main__":
    n = 15
    A, b, h = assemble_poisson(n)
    uh = spsolve(A, b)
    residual = np.linalg.norm(b - A @ uh, np.inf)
    mid = n // 2  # x = y = 0.5 is grid index (n-1)/2 for odd n
    centre = uh[mid * n + mid]
    print(f"n = {n}: algebraic residual ||b - A u_h||_inf = {residual:.3e}")
    print(f"centre u_h(0.5, 0.5) = {centre:.4f}  (expect ~1.003, u* = 1)")
    assert residual < 1e-10, "residual not near machine precision"
    assert abs(centre - 1.003) < 2e-3, "centre value off the P1 sanity value"
    print("environment OK: sparse assembly + spsolve work")
