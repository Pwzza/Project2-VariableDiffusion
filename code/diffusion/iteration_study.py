"""Expanded unpreconditioned CG sweep with independently checked residuals."""
import csv
import json
from time import perf_counter
import numpy as np
from scipy.sparse.linalg import cg
from .models import Model
from .solvers import assemble_fdm,assemble_fem


def run(root,quick=False):
    ns=(7,15,31) if quick else (15,31,63,127)
    contrasts=(1.,10.,100.) if quick else (1.,3.,10.,30.,100.)
    rows=[]
    for rho in contrasts:
        model=Model('layered',rho)
        for N in ns:
            for method in ('fdm_resistance','fem_fitted'):
                s=assemble_fdm(model,N,'resistance') if method=='fdm_resistance' else assemble_fem(model,N)
                count=[0]
                def callback(x):count[0]+=1
                start=perf_counter()
                u,info=cg(s.A,s.b,x0=np.zeros_like(s.b),rtol=9e-11,atol=0.,maxiter=20000,callback=callback)
                seconds=perf_counter()-start
                residual=float(np.linalg.norm(s.b-s.A@u)/np.linalg.norm(s.b))
                if info!=0 or residual>1e-10:
                    raise AssertionError(f'{method}, N={N}, rho={rho}: info={info}, residual={residual}')
                coords=np.meshgrid(s.x,s.y,indexing='ij') if s.basis is None else s.basis.doflocs
                exact=model.exact(*coords).ravel()[s.free]
                rows.append(dict(method=method,N=N,h=1/N,contrast=rho,unknowns=len(s.b),
                                 iterations=count[0],solver_info=int(info),rtol=9e-11,
                                 acceptance_relative_residual=1e-10,residual_rel2=residual,
                                 nodal_linf=float(np.max(abs(u-exact))),seconds=seconds,
                                 preconditioner='none',initial_guess='zero'))
                print(f'CG {method} N={N} contrast={rho:g}: {count[0]} iterations; residual={residual:.2e}',flush=True)
    (root/'data'/'iterations.json').write_text(json.dumps(rows,indent=2)+'\n',encoding='utf-8')
    with (root/'data'/'iterations.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return rows
