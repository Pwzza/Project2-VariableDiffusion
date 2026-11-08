import json
import numpy as np
from .models import Model
from .solvers import assemble_fdm,assemble_fem,solve

STYLES={'fdm_arithmetic':dict(color='#d07919',marker='o',linestyle='-'),
        'fdm_resistance':dict(color='#2176ae',marker='s',linestyle='--',markerfacecolor='none'),
        'fem_fitted':dict(color='#8054a0',marker='x',linestyle=':')}
LABELS={'fdm_arithmetic':'FDM arithmetic','fdm_resistance':'FDM resistance','fem_fitted':'FEM fitted'}

def benchmark(root):
    rows=[]
    for model_name in ('smooth','layered'):
        cases=json.loads((root/'data'/f'{model_name}.json').read_text())
        for case in cases:
            def create():
                m=Model(model_name);n=case['N']
                if case['method']=='fem_fitted': return assemble_fem(m,n)
                return assemble_fdm(m,n,case['method'].removeprefix('fdm_'))
            solve(create())
            times=[]
            for _ in range(3):
                s=solve(create());times.append(s.assembly_s+s.solve_s)
            rows.append({**{k:case[k] for k in ('model','method','N','unknowns','l2','nodal_linf')},
                         'seconds':times,'median_seconds':float(np.median(times)),
                         'flux_relative_linf':case.get('flux_relative_linf')})
            print(f"cost: {model_name} {case['method']} N={case['N']} {np.median(times):.3g}s",flush=True)
    (root/'data'/'cost_accuracy.json').write_text(json.dumps(rows,indent=2)+'\n')

def figures(root):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from scipy.interpolate import RegularGridInterpolator
    from skfem import MeshTri,Basis,ElementTriP1
    def read(name):return json.loads((root/'data'/f'{name}.json').read_text())
    def save(fig,name):
        fig.tight_layout()
        for ext in ('png','pdf'):fig.savefig(root/'figures'/f'{name}.{ext}',bbox_inches='tight',dpi=190)
        plt.close(fig)
    rows=read('layered');N=min(r['N'] for r in rows)
    fig,axs=plt.subplots(1,2,figsize=(11,4))
    for k,(method,style) in enumerate(STYLES.items()):
        p=read(f'layered_{method}_N{N}_profile');x=np.array(p['temperature_x'])
        for ax,xx,yy in ((axs[0],x,np.array(p['temperature'])-Model('layered').exact(x,0)),
                         (axs[1],p['x'],p['q'])):
            ax.plot(xx,yy,**style,markevery=(k*2,max(4,len(xx)//12)),label=LABELS[method],lw=1.8)
    axs[1].axhline(-20/11,color='.4',ls='-.',lw=.8,label='Exact flux')
    for ax,title in zip(axs,('Temperature error','Signed flux $q_x$')):
        ax.axvline(.5,color='.7',ls=':');ax.set(title=title,xlabel='x');ax.legend(fontsize=9);ax.grid(alpha=.2)
    save(fig,'slides_layered_result')
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for method,style in STYLES.items():
        rs=[r for r in rows if r['method']==method]
        for ax,key in zip(axs,('l2','nodal_linf')):
            ax.loglog([r['h'] for r in rs],[r[key] for r in rs],**style,label=LABELS[method])
            ax.set(xlabel='$h=1/N$',ylabel='Continuous $L^2$ error' if key=='l2' else 'Nodal maximum error')
            ax.legend(fontsize=9);ax.grid(alpha=.2)
    for ax in axs:
        hs=[r['h'] for r in rows if r['method']=='fdm_arithmetic']
        ax.set_xticks(hs,labels=[f'1/{round(1/h)}' for h in hs])
        ax.xaxis.set_minor_formatter(NullFormatter())
    save(fig,'layered_convergence')

    import matplotlib.tri as mtri
    error_fields=[]
    for method in STYLES:
        d=np.load(root/'data'/f'layered_{method}_N{N}.npz')
        if method.startswith('fdm'):
            X,Y=np.meshgrid(d['x'],d['y'],indexing='ij')
            error=d['u']-Model('layered').exact(X,Y)
        else:
            error=d['u']-Model('layered').exact(*d['points'])
        error_fields.append((method,d,error))
    limit=max(np.max(abs(error)) for _,_,error in error_fields)
    fig,axs=plt.subplots(1,3,figsize=(12,3.8),layout='constrained')
    for ax,(method,d,error) in zip(axs,error_fields):
        if method.startswith('fdm'):
            im=ax.pcolormesh(d['x'],d['y'],error.T,shading='auto',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
        else:
            tri=mtri.Triangulation(*d['points'],d['triangles'].T)
            im=ax.tripcolor(tri,error,shading='gouraud',cmap='RdBu_r',vmin=-limit,vmax=limit,rasterized=True)
        p=read(f'layered_{method}_N{N}_profile')
        ax.axvline(.5,ls=':',color='.4',lw=1)
        ax.axhline(p['slice_y'],ls='--',color='.25',lw=1)
        ax.set(title=LABELS[method],xlabel='x',ylabel='y',aspect='equal',xlim=(0,1),ylim=(0,1))
    fig.colorbar(im,ax=axs,label='$u_h-u^*$',shrink=.8)
    for ext in ('png','pdf'):fig.savefig(root/'figures'/f'layered_error_maps.{ext}',bbox_inches='tight',dpi=190)
    plt.close(fig)
    # Actual numerical q fields, FDM interpolates shared face fluxes
    # FEM samples the elementwise gradient, the common scale preserves magnitude
    xx,yy=np.meshgrid(np.linspace(.06,.94,10),np.linspace(.06,.94,10),indexing='ij')
    points=np.array([xx.ravel(),yy.ravel()]);fig,axs=plt.subplots(1,3,figsize=(12,4))
    vector_stats=[]
    for ax,(method,style) in zip(axs,STYLES.items()):
        d=np.load(root/'data'/f'layered_{method}_N{N}.npz')
        if method.startswith('fdm'):
            s=assemble_fdm(Model('layered'),N,method.removeprefix('fdm_'))
            qx=-s.ax*np.diff(d['u'],axis=0)*N;qy=-s.ay*np.diff(d['u'],axis=1)*N
            q=np.array([RegularGridInterpolator(((s.x[:-1]+s.x[1:])/2,s.y),qx)(points.T),
                        RegularGridInterpolator((s.x,(s.y[:-1]+s.y[1:])/2),qy)(points.T)])
        else:
            mesh=MeshTri(d['points'],d['triangles']);basis=Basis(mesh,ElementTriP1())
            cells=mesh.element_finder()(*points)
            q=-Model('layered').a(*points)*basis.interpolate(d['u']).grad[:,cells,0]
        magnitude=np.linalg.norm(q,axis=0)
        vector_stats.append(dict(method=method,N=N,min_magnitude=float(magnitude.min()),
                                 max_magnitude=float(magnitude.max()),max_abs_qy=float(abs(q[1]).max())))
        ax.quiver(*points,*q,angles='xy',scale_units='xy',scale=25,pivot='middle',color=style['color'],width=.004)
        ax.axvline(.5,ls=':',color='.5');ax.set(xlim=(0,1),ylim=(0,1),aspect='equal',xlabel='x',ylabel='y',title=LABELS[method])
    save(fig,'layered_numerical_flux_vectors')
    (root/'data'/'flux_vector_inspection.json').write_text(json.dumps(vector_stats,indent=2)+'\n')
    if not (root/'data'/'cost_accuracy.json').exists():return
    cost=read('cost_accuracy')
    for model in ('smooth','layered'):
        fig,axs=plt.subplots(1,2,figsize=(10,4))
        for method,style in STYLES.items():
            rs=[r for r in cost if r['model']==model and r['method']==method]
            if not rs:continue
            for ax,key in zip(axs,('l2','nodal_linf' if model=='smooth' else 'flux_relative_linf')):
                ax.loglog([r['median_seconds'] for r in rs],[r[key] for r in rs],**style,label=LABELS[method])
                ax.set(xlabel='Median assembly + direct solve (s)',ylabel={'l2':'Continuous $L^2$ error','nodal_linf':'Nodal maximum error','flux_relative_linf':'Relative slice flux error'}[key]);ax.legend(fontsize=9);ax.grid(alpha=.2)
        save(fig,f'{model}_cost_accuracy')
