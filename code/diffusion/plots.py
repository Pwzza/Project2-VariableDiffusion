"""Headless presentation figures, regenerated exclusively from saved results.

Analytic introduction/method diagrams are mathematical illustrations; numerical
result plots read data/*.json and *.npz. All assets export PDF and PNG.
"""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.ticker import NullFormatter
from .models import Model

COLORS={'fdm_arithmetic':'#d07919','fdm_resistance':'#2176ae','fem_fitted':'#8054a0'}
LABELS={'fdm_arithmetic':'FDM arithmetic','fdm_resistance':'FDM resistance','fem_fitted':'FEM fitted'}


def generate(root):
    plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,
                         'savefig.dpi':190,'pdf.fonttype':42})
    def read(name): return json.loads((root/'data'/f'{name}.json').read_text())
    def save(fig,name):
        fig.tight_layout()
        for ext in ('png','pdf'): fig.savefig(root/'figures'/f'{name}.{ext}',bbox_inches='tight')
        plt.close(fig)
    x=np.linspace(0,1,401);X,Y=np.meshgrid(x,x,indexing='ij')
    fig,axs=plt.subplots(2,2,figsize=(9,6))
    for j,m in enumerate((Model(),Model('layered'))):
        for i,(field,title) in enumerate(((m.a(X,Y),'Conductivity'),(m.exact(X,Y),'Exact temperature'))):
            im=axs[i,j].pcolormesh(X,Y,field,shading='auto',cmap='viridis',rasterized=True)
            axs[i,j].set(xlabel='x',ylabel='y',aspect='equal',title=f'{m.name}: {title}')
            fig.colorbar(im,ax=axs[i,j])
    save(fig,'intro_models')
    fig,axs=plt.subplots(3,1,figsize=(8,6),sharex=True)
    m=Model('layered')
    for ax,v,title in zip(axs,[m.a(x,0),m.exact(x,0),m.flux(x,0)[0]],
                          ['Conductivity $a$','Temperature $u^*$','Flux $q_x^*$']):
        ax.plot(x,v,lw=2);ax.axvline(.5,color='gray',ls=':');ax.set_ylabel(title)
    axs[-1].set_xlabel('x');save(fig,'intro_interface')
    # Explicit methodology diagrams requested for the presentation.
    fig,axs=plt.subplots(1,2,figsize=(9,4))
    ax=axs[0]
    for t in (-1,0,1):
        ax.plot([-1,1],[t,t],color='.8');ax.plot([t,t],[-1,1],color='.8')
    ax.scatter([0,1,-1,0,0],[0,0,0,1,-1],color=COLORS['fdm_resistance'])
    for dx,dy,label in ((.7,0,'E'),(-.7,0,'W'),(0,.7,'N'),(0,-.7,'S')):
        ax.annotate('',xy=(dx,dy),xytext=(0,0),arrowprops={'arrowstyle':'->','lw':2})
        ax.text(dx+.06,dy+.06,label)
    ax.set(title='FDM: shared face fluxes',aspect='equal');ax.axis('off')
    from skfem import MeshTri
    mesh=MeshTri.init_tensor(np.array([0,.25,.5,.75,1]),np.linspace(0,1,5))
    axs[1].triplot(mesh.p[0],mesh.p[1],mesh.t.T,color=COLORS['fem_fitted'])
    axs[1].axvline(.5,lw=3,color='black');axs[1].set(title='FEM: interface-fitted triangles',aspect='equal',xlabel='x',ylabel='y')
    save(fig,'methods_geometry')
    fig=plt.figure(figsize=(6,4));ax=fig.add_subplot(111,projection='3d')
    phi=np.zeros(mesh.p.shape[1]);phi[np.argmin(np.sum((mesh.p-.5)**2,axis=0))]=1.
    ax.plot_trisurf(mesh.p[0],mesh.p[1],phi,triangles=mesh.t.T,color=COLORS['fem_fitted'],
                    edgecolor='white',linewidth=.5,alpha=.9)
    ax.set(xlabel='x',ylabel='y',zlabel=r'$\phi_i$',title='P1 nodal basis: linear on each triangle')
    save(fig,'fem_basis')
    fig,ax=plt.subplots(figsize=(8,3))
    ax.plot([0,.5,1],[0,0,0],lw=8,color='.8')
    ax.scatter([0,.5,1],[0,0,0],s=80,color='black',zorder=3)
    ax.text(.25,.12,'$d_L$, $a_L$',ha='center');ax.text(.75,.12,'$d_R$, $a_R$',ha='center')
    ax.text(.5,-.25,r'$R=d_L/a_L+d_R/a_R,\qquad a_{\mathrm{edge}}=h/R$',ha='center',fontsize=16)
    ax.annotate('constant flux',xy=(.8,.4),xytext=(.2,.4),arrowprops={'arrowstyle':'->'},va='center')
    ax.set(xlim=(-.05,1.05),ylim=(-.45,.65),title='Interface edge: thermal resistances in series');ax.axis('off')
    save(fig,'interface_resistance')
    from scipy.sparse import load_npz
    fig,axs=plt.subplots(1,2,figsize=(8,4))
    for ax,method in zip(axs,('fdm_arithmetic','fem_fitted')):
        A=load_npz(root/'data'/f'smooth_{method}_matrix.npz')
        ax.spy(A,markersize=1,color=COLORS[method]);ax.set_title(LABELS[method])
    save(fig,'matrix_structure')
    smooth=read('smooth');fits=read('fits')
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for method in ('fdm_arithmetic','fem_fitted'):
        rows=[r for r in smooth if r['method']==method]
        h=np.array([r['h'] for r in rows])
        for ax,norm,label in zip(axs,('l2','nodal_linf'),('Continuous $L^2$ error','Nodal maximum error')):
            fit=next(f for f in fits if f['method']==method and f['norm']==norm)
            ax.loglog(h,[r[norm] for r in rows],'o-',color=COLORS[method],label=f"{LABELS[method]}, p={fit['order']:.2f}")
            ax.set(xlabel='$h=1/N$',ylabel=label);ax.grid(alpha=.2,which='both')
    h=np.array([r['h'] for r in smooth if r['method']=='fdm_arithmetic'])
    for ax,norm in zip(axs,('l2','nodal_linf')):
        e=next(r[norm] for r in smooth)
        ax.loglog(h,e*(h/h[0])**2,'k--',alpha=.5,label='$O(h^2)$ reference');ax.legend(fontsize=9)
        ax.set_xticks(h,labels=[f'1/{round(1/v)}' for v in h])
        ax.xaxis.set_minor_formatter(NullFormatter())
    save(fig,'smooth_convergence')
    # Same color limits for error maps, actual nodal errors.
    N=max(r['N'] for r in smooth); fig,axs=plt.subplots(1,2,figsize=(9,4))
    cases=[]
    for method in ('fdm_arithmetic','fem_fitted'):
        d=np.load(root/'data'/f'smooth_{method}_N{N}.npz')
        if method.startswith('fdm'):
            xx,yy=np.meshgrid(d['x'],d['y'],indexing='ij');err=d['u']-Model().exact(xx,yy)
        else: err=d['u']-Model().exact(*d['points'])
        cases.append((method,d,err))
    lim=max(np.max(abs(e)) for _,_,e in cases)
    for ax,(method,d,e) in zip(axs,cases):
        if method.startswith('fdm'): im=ax.pcolormesh(d['x'],d['y'],e.T,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='auto')
        else: im=ax.tripcolor(mtri.Triangulation(*d['points'],d['triangles'].T),e,cmap='RdBu_r',vmin=-lim,vmax=lim,shading='gouraud')
        ax.set(title=LABELS[method],xlabel='x',ylabel='y',aspect='equal');fig.colorbar(im,ax=ax)
    save(fig,'smooth_error_maps')
    layered=read('layered');N=min(r['N'] for r in layered)
    fig,axs=plt.subplots(1,3,figsize=(13,4))
    for method in COLORS:
        p=read(f'layered_{method}_N{N}_profile')
        tx=np.array(p['temperature_x']);u=np.array(p['temperature'])
        label=f"{LABELS[method]} (y={p['slice_y']:.3f})"
        axs[0].plot(tx,u,'.-',color=COLORS[method],label=label)
        axs[1].plot(p['x'],p['q'],'.-',color=COLORS[method])
        axs[2].plot(tx,u-Model('layered').exact(tx,0),'.-',color=COLORS[method])
    axs[0].plot(x,m.exact(x,0),'k--',label='Exact');axs[0].legend(fontsize=8)
    axs[1].axhline(-20/11,color='black',ls='--')
    for ax,title in zip(axs,('Temperature','Signed normal flux $q_x$','Temperature error on slice')):
        ax.axvline(.5,color='gray',ls=':');ax.set(xlabel='x',title=title);ax.grid(alpha=.2)
    save(fig,'layered_profiles')
    fig,axs=plt.subplots(1,2,figsize=(10,4))
    for method in COLORS:
        rows=[r for r in layered if r['method']==method]
        for ax,key,title in zip(axs,('nodal_linf','flux_relative_linf'),('Nodal maximum error','Relative flux error')):
            ax.loglog([r['h'] for r in rows],[max(r[key],1e-16) for r in rows],'o-',color=COLORS[method],label=LABELS[method])
            ax.set(xlabel='$h=1/N$',ylabel=title);ax.grid(which='both',alpha=.2);ax.legend(fontsize=9)
    h=np.array([r['h'] for r in layered if r['method']=='fdm_arithmetic'])
    for ax in axs:
        ax.set_xticks(h,labels=[f'1/{round(1/v)}' for v in h])
        ax.xaxis.set_minor_formatter(NullFormatter())
    save(fig,'layered_errors')
    fig,axs=plt.subplots(1,2,figsize=(9,4))
    iterations=read('iterations')
    for ax,method in zip(axs,('fdm_resistance','fem_fitted')):
        for rho in (1.,10.,100.):
            rows=[r for r in iterations if r['method']==method and r['contrast']==rho]
            ax.plot([r['N'] for r in rows],[r['iterations'] for r in rows],'o-',label=f'contrast {rho:g}')
        ax.set(title=LABELS[method],xlabel='Intervals N',ylabel='Unpreconditioned CG iterations');ax.legend(fontsize=9);ax.grid(alpha=.2)
    save(fig,'contrast_iterations')
