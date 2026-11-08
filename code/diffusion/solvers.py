from dataclasses import dataclass
from time import perf_counter
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from skfem import Basis, MeshTri, ElementTriP1, BilinearForm, LinearForm, asm
from skfem.helpers import dot, grad


@dataclass
class Solution:
    model: object
    method: str
    N: int
    A: object
    b: np.ndarray
    u: np.ndarray
    free: np.ndarray
    x: np.ndarray
    y: np.ndarray
    assembly_s: float
    solve_s: float = 0.
    basis: object = None
    ax: object = None
    ay: object = None

    def with_free(self, values):
        u = self.u.copy()
        u.ravel()[self.free] = values
        return u


def face_coefficient(model, x0, x1, y0, y1, averaging):
    a0, a1 = model.a(x0, y0), model.a(x1, y1)
    if averaging == 'arithmetic':
        return (a0+a1)/2
    if averaging != 'resistance':
        raise ValueError('Unknown face averaging')
    if model.name == 'smooth':
        return (a0+a1)/2
    # Exact resistance fraction also handles an interface at a node
    if x1 == x0:
        return model.a(x0, (y0+y1)/2)
    left_fraction = np.clip((.5-x0)/(x1-x0), 0., 1.)
    return 1./(left_fraction+(1-left_fraction)/model.contrast)


def assemble_fdm(model, N, averaging='arithmetic'):
    if N < 3:
        raise ValueError('N must be at least 3 intervals')
    start = perf_counter()
    x = y = np.linspace(0, 1, N+1)
    X, Y = np.meshgrid(x, y, indexing='ij')
    ax = np.empty((N, N+1)); ay = np.empty((N+1, N))
    for i in range(N):
        ax[i] = face_coefficient(model, x[i], x[i+1], y, y, averaging)
    for i in range(N+1):
        ay[i] = face_coefficient(model, x[i], x[i], y[:-1], y[1:], averaging)
    u = model.exact(X, Y).copy()
    u[1:-1, 1:-1] = 0.
    ids = -np.ones_like(u, dtype=int)
    ids[1:-1, 1:-1] = np.arange((N-1)**2).reshape(N-1, N-1)
    b = model.forcing(X[1:-1, 1:-1], Y[1:-1, 1:-1]).ravel().copy()
    rows, cols, vals = [], [], []
    for i in range(1, N):
        for j in range(1, N):
            k = ids[i,j]; diagonal = 0.
            for ii,jj,a in ((i+1,j,ax[i,j]), (i-1,j,ax[i-1,j]),
                            (i,j+1,ay[i,j]), (i,j-1,ay[i,j-1])):
                w = float(a)*N*N
                diagonal += w
                if ids[ii,jj] >= 0:
                    rows.append(k); cols.append(ids[ii,jj]); vals.append(-w)
                else:
                    b[k] += w*u[ii,jj]
            rows.append(k); cols.append(k); vals.append(diagonal)
    A = coo_matrix((vals,(rows,cols)), shape=((N-1)**2,)*2).tocsr()
    free = np.flatnonzero(ids.ravel() >= 0)
    return Solution(model, 'fdm_'+averaging, N, A, b, u, free, x, y,
                    perf_counter()-start, ax=ax, ay=ay)


def assemble_fem(model, N, quadrature=6):
    start = perf_counter()
    x = np.linspace(0, 1, N+1)
    if model.name == 'layered':
        x = np.unique(np.r_[x, .5])
    y = np.linspace(0, 1, N+1)
    mesh = MeshTri.init_tensor(x, y)
    basis = Basis(mesh, ElementTriP1(), intorder=quadrature)

    @BilinearForm
    def stiffness(u, v, w):
        return model.a(w.x[0], w.x[1])*dot(grad(u), grad(v))

    @LinearForm
    def load(v, w):
        return model.forcing(w.x[0], w.x[1])*v

    K = asm(stiffness, basis).tocsr()
    F = asm(load, basis)
    boundary = basis.get_dofs().all()
    free = np.setdiff1d(np.arange(basis.N), boundary)
    u = np.zeros(basis.N)
    u[boundary] = model.exact(*basis.doflocs[:,boundary])
    b = F[free]-K[free][:,boundary]@u[boundary]
    return Solution(model, 'fem_fitted', N, K[free][:,free].tocsr(), b,
                    u, free, x, y, perf_counter()-start, basis=basis)


def solve(system):
    start = perf_counter()
    system.u = system.with_free(spsolve(system.A, system.b))
    system.solve_s = perf_counter()-start
    return system
