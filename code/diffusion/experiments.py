"""Reproducible experiments; machine-readable data precedes plotting."""
import csv
import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from time import perf_counter
import numpy as np
from scipy.sparse import save_npz
from scipy.sparse.linalg import cg, eigsh
from .models import Model
from .solvers import assemble_fdm, assemble_fem, solve
from .diagnostics import metrics, physical_errors, flux_profile


def write_json(path, obj):
    path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')


def save_case(root,s,row):
    stem=f'{s.model.name}_{s.method}_N{s.N}'
    data=dict(u=s.u,x=s.x,y=s.y)
    if s.basis is not None:
        data.update(points=s.basis.mesh.p,triangles=s.basis.mesh.t)
    np.savez_compressed(root/'data'/f'{stem}.npz',**data)
    if s.model.name=='layered':
        p=flux_profile(s)
        write_json(root/'data'/f'{stem}_profile.json',
                   {k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in p.items()})


def add_orders(rows):
    fits=[]
    for method in sorted(set(r['method'] for r in rows)):
        group=[r for r in rows if r['method']==method]
        for norm in ('l2','nodal_linf','energy'):
            valid=[r for r in group if norm in r]
            if len(valid)<2: continue
            # Exclude roundoff-level errors from fitting; do not invent a slope.
            for prev,current in zip(valid[:-1],valid[1:]):
                if min(prev[norm],current[norm])>1e-11:
                    current[norm+'_pair_order']=float(np.log(prev[norm]/current[norm])/np.log(prev['h']/current['h']))
            window=valid[-3:]
            if min(r[norm] for r in window)<=1e-11: continue
            x=np.log([r['h'] for r in window]); y=np.log([r[norm] for r in window])
            c=np.polyfit(x,y,1)
            fits.append(dict(method=method,norm=norm,order=float(c[0]),
                             fit_N=[r['N'] for r in window],
                             max_log_deviation=float(np.max(np.abs(y-np.polyval(c,x))))))
    return fits


def independent_forcing_check():
    """Differentiate analytic FLUX numerically, not the supplied forcing formula."""
    m=Model(); x=np.linspace(.08,.92,13); y=np.linspace(.09,.91,13)
    x,y=np.meshgrid(x,y,indexing='ij'); d=1e-5
    divergence=((m.flux(x+d,y)[0]-m.flux(x-d,y)[0])+
                (m.flux(x,y+d)[1]-m.flux(x,y-d)[1]))/(2*d)
    defect=float(np.max(np.abs(divergence-m.forcing(x,y))))
    if defect>1e-7: raise AssertionError(f'Forcing check failed: {defect}')
    return dict(method='centered numerical divergence of analytic flux',step=d,
                maximum_absolute_defect=defect,threshold=1e-7)


def run(root, quick=False):
    for folder in ('data','figures','notes/evidence'):
        (root/folder).mkdir(parents=True,exist_ok=True)
    start=perf_counter()
    smooth_N=[8,16,32,64] if quick else [16,32,64,128]
    layered_N=[7,15,31,63] if quick else [31,63,127,255]
    checks={'forcing':independent_forcing_check()}
    smooth=[]; layered=[]
    for model,Ns,rows in ((Model(),smooth_N,smooth),(Model('layered'),layered_N,layered)):
        for N in Ns:
            constructors=[lambda:assemble_fdm(model,N,'arithmetic'),lambda:assemble_fem(model,N)]
            if model.name=='layered':
                constructors.insert(1,lambda:assemble_fdm(model,N,'resistance'))
            for constructor in constructors:
                s=solve(constructor()); row=metrics(s); rows.append(row)
                save_case(root,s,row)
                if N==Ns[0]:
                    save_npz(root/'data'/f'{model.name}_{s.method}_matrix.npz',s.A)
                print(f"{model.name:7} {s.method:15} N={N:3} L2={row['l2']:.4e} "
                      f"max={row['nodal_linf']:.3e}",flush=True)
    # Quantify quadrature sensitivity and residual correction on finest smooth mesh.
    N=smooth_N[-1]
    base=solve(assemble_fem(Model(),N,6)); high=solve(assemble_fem(Model(),N,8))
    delta=physical_errors(base,high.u-base.u,exact=False)['l2']
    fem_error=[r['l2'] for r in smooth if r['method']=='fem_fitted' and r['N']==N][0]
    checks['quadrature']=dict(N=N,orders=[6,8],solution_difference_l2=delta,
                              difference_to_discretization_ratio=delta/fem_error)
    # Positive-definiteness numerical corroboration on a small representative mesh.
    checks['smallest_eigenvalues']=[]
    for model in (Model(),Model('layered')):
        for assemble in (lambda:assemble_fdm(model,15,'resistance'),lambda:assemble_fem(model,15)):
            s=assemble()
            eigenvalue=float(eigsh(s.A,k=1,which='SA',return_eigenvectors=False,tol=1e-9)[0])
            checks['smallest_eigenvalues'].append(dict(model=model.name,method=s.method,N=15,value=eigenvalue))
    iterations=[]
    for rho in (1.,10.,100.):
        for N in (7,15,31) if quick else (15,31,63):
            for assemble in (lambda:assemble_fdm(Model('layered',rho),N,'resistance'),
                             lambda:assemble_fem(Model('layered',rho),N)):
                s=assemble(); counter=[0]
                def callback(x): counter[0]+=1
                t=perf_counter()
                u,info=cg(s.A,s.b,rtol=1e-10,atol=0.,maxiter=20000,callback=callback)
                elapsed=perf_counter()-t
                residual=float(np.linalg.norm(s.b-s.A@u)/np.linalg.norm(s.b))
                if info!=0 or residual>1.1e-10:
                    raise AssertionError(f'CG did not attain target: info={info}, residual={residual}')
                iterations.append(dict(method=s.method,N=N,h=1/N,contrast=rho,
                                       iterations=counter[0],rtol=1e-10,residual_rel2=residual,
                                       seconds=elapsed,preconditioner='none',initial_guess='zero'))
    fits=add_orders(smooth)
    checks['acceptance']={
        'smooth_l2_orders':all(1.8<f['order']<2.2 for f in fits if f['norm']=='l2'),
        'smooth_algebraic_separation':all(r['algebraic_l2_correction']<1e-3*r['l2'] for r in smooth),
        'quadrature_separation':delta<1e-3*fem_error,
        'symmetric':all(r['symmetry_rel']<1e-12 for r in smooth+layered),
        'positive_eigenvalues':all(r['value']>0 for r in checks['smallest_eigenvalues']),
        'layered_exact_nodes':all(r['nodal_linf']<1e-9 for r in layered if r['method']!='fdm_arithmetic'),
        'layered_exact_flux':all(r['flux_relative_linf']<1e-8 for r in layered if r['method']!='fdm_arithmetic'),
        'layered_bounds':all(r['boundary_bound_violation']<1e-10 for r in layered)}
    for name,rows in (('smooth',smooth),('layered',layered),('iterations',iterations)):
        write_json(root/'data'/f'{name}.json',rows)
        keys=sorted(set().union(*(r.keys() for r in rows)))
        with (root/'data'/f'{name}.csv').open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=keys);writer.writeheader();writer.writerows(rows)
    write_json(root/'data'/'fits.json',fits)
    write_json(root/'data'/'checks.json',checks)
    metadata=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),quick=quick,
                  smooth_N=smooth_N,layered_N=layered_N,python=sys.version,platform=platform.platform(),
                  packages={p:version(p) for p in ('numpy','scipy','matplotlib','scikit-fem')},
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted((root/'code'/'diffusion').glob('*.py'))},
                  timing_note='Single-run wall times; diagnostic only, not a controlled benchmark.',
                  norm_note='FDM L2 uses interface-aware resistance reconstruction for layered model; l2_bilinear is separate. FEM L2 is integrated P1. Nodal Linf is not continuous Linf.',
                  total_experiment_seconds=perf_counter()-start)
    write_json(root/'data'/'metadata.json',metadata)
    lines=['# Generated numerical evidence','',
           f"Run mode: {'quick (development)' if quick else 'full'}. See data/metadata.json.",'',
           '## Smooth model: selected fits','',
           '| Method | Norm | Order | Fitted N | Max log deviation |',
           '|---|---|---:|---|---:|']
    for f in fits:
        lines.append(f"| {f['method']} | {f['norm']} | {f['order']:.5f} | {f['fit_N']} | {f['max_log_deviation']:.3g} |")
    lines += ['','## Layered model: finest mesh','',
              '| Method | N | Nodal max error | Relative slice flux error | Boundary violation |',
              '|---|---:|---:|---:|---:|']
    for r in layered:
        if r['N']==layered_N[-1]:
            lines.append(f"| {r['method']} | {r['N']} | {r['nodal_linf']:.5g} | {r['flux_relative_linf']:.5g} | {r['boundary_bound_violation']:.3g} |")
    lines += ['','## Acceptance gates','']
    lines += [f"- {name}: {'PASS' if passed else 'FAIL'}" for name,passed in checks['acceptance'].items()]
    lines += ['','Residual correction is a numerical algebraic-error estimate, not a certified bound.',
              'Timings are single-run diagnostics; no controlled speed ranking is claimed.','']
    (root/'notes'/'evidence'/'numerical_results.md').write_text('\n'.join(lines),encoding='utf-8')
    if not all(checks['acceptance'].values()):
        raise AssertionError(f"Acceptance failure; inspect data/checks.json: {checks['acceptance']}")
    return checks
