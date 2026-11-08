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

    qgrid=np.linspace(.04,.96,15)
    QX,QY=np.meshgrid(qgrid,qgrid,indexing='ij')
    models=(Model(),Model('layered'))
    fig,axs=plt.subplots(1,2,figsize=(9,4.2),layout='constrained')
    for ax,model in zip(axs,models):
        arrow_scale=np.max(np.linalg.norm(model.flux(QX,QY),axis=0))/.05
        im=ax.pcolormesh(X,Y,model.exact(X,Y),shading='auto',cmap='Blues',
                         vmin=0,vmax=1,alpha=.35,rasterized=True)
        q=model.flux(QX,QY)
        arrows=ax.quiver(QX,QY,q[0],q[1],angles='xy',scale_units='xy',
                         scale=arrow_scale,pivot='middle',width=.004,
                         headwidth=3.5,headlength=4.5,color='#182b49')
        ax.quiverkey(arrows,.78,1.06,2.,r'$|\mathbf{q}|=2$',labelpos='E',coordinates='axes')
        if model.name=='layered': ax.axvline(.5,color='#b05a12',ls='--',lw=1.5)
        ax.set(xlabel='x',ylabel='y',title=model.name.capitalize(),aspect='equal',xlim=(0,1),ylim=(0,1))
    fig.colorbar(im,ax=axs,label='Exact temperature',shrink=.8)

    for ext in ('png','pdf'):
        fig.savefig(root/'figures'/f'flux_vector_field.{ext}',bbox_inches='tight')
    plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(8,3.2))
    for ax,model in zip(axs,models):
        im=ax.pcolormesh(X,Y,model.a(X,Y),shading='auto',cmap='viridis',rasterized=True)
        ax.set(title=model.name.capitalize(),xlabel='x',ylabel='y',aspect='equal')
        fig.colorbar(im,ax=ax,label='$a$')
    save(fig,'slides_conductivity')
    fig,axs=plt.subplots(3,1,figsize=(8,6),sharex=True)
    m=Model('layered')
    for ax,v,title in zip(axs,[m.a(x,0),m.exact(x,0),m.flux(x,0)[0]],
                          ['Conductivity $a$','Temperature $u^*$','Flux $q_x^*$']):
        ax.plot(x,v,lw=2);ax.axvline(.5,color='gray',ls=':');ax.set_ylabel(title)
    axs[-1].set_xlabel('x');save(fig,'intro_interface')

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
    axs[1].remove()
    ax=fig.add_subplot(122,projection='3d')
    phi=np.zeros(mesh.p.shape[1]);phi[np.argmin(np.sum((mesh.p-.5)**2,axis=0))]=1.
    ax.plot_trisurf(mesh.p[0],mesh.p[1],phi,triangles=mesh.t.T,
                    color=COLORS['fem_fitted'],edgecolor='white',linewidth=.5,alpha=.9)
    ax.set(xlabel='x',ylabel='y',zlabel=r'$\phi_i$',title='FEM: P1 nodal basis')
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
    fig,ax=plt.subplots(figsize=(6.2,3.7))
    for method in ('fdm_arithmetic','fem_fitted'):
        rows=[r for r in smooth if r['method']==method]
        h=np.array([r['h'] for r in rows]);e=np.array([r['l2'] for r in rows])
        p=next(f['order'] for f in fits if f['method']==method and f['norm']=='l2')
        ax.loglog(h,e,'o-',color=COLORS[method],label=f"{LABELS[method]}, p={p:.3f}")
    ref=np.array([r['l2'] for r in smooth if r['method']=='fdm_arithmetic'])[0]
    ax.loglog(h,ref*(h/h[0])**2,'k--',alpha=.5,label='$O(h^2)$ reference')
    ax.set(xlabel='$h=1/N$',ylabel='Continuous $L^2$ error')
    ax.set_xticks(h,labels=[f'1/{round(1/v)}' for v in h]);ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(fontsize=10);ax.grid(which='both',alpha=.2)
    save(fig,'slides_smooth_l2')

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
    fig,axs=plt.subplots(1,2,figsize=(9,3.5))
    for method in COLORS:
        p=read(f'layered_{method}_N{N}_profile')
        tx=np.array(p['temperature_x']);u=np.array(p['temperature'])
        axs[0].plot(tx,u-Model('layered').exact(tx,0),'.-',color=COLORS[method],label=LABELS[method])
        axs[1].plot(p['x'],p['q'],'.-',color=COLORS[method])
    axs[1].axhline(-20/11,color='black',ls='--',label='Exact flux');axs[1].legend(fontsize=9)
    axs[0].legend(fontsize=9)
    for ax,title in zip(axs,('Temperature error','Signed flux $q_x$')):
        ax.axvline(.5,color='gray',ls=':');ax.set(xlabel='x',title=title);ax.grid(alpha=.2)
    save(fig,'slides_layered_result')

    tex=[]
    for method,command in (('fdm_arithmetic','FDMOrder'),('fem_fitted','FEMOrder')):
        p=next(f['order'] for f in fits if f['method']==method and f['norm']=='l2')
        tex.append(f'{chr(92)}newcommand{{{chr(92)}{command}}}{{{p:.4f}}}')
    iterations=read('iterations')
    for rho in (1,10,100):
        for method,command in (('fdm_resistance','FD'),('fem_fitted','FE')):
            r=next(r for r in iterations if r['contrast']==rho and r['method']==method and r['N']==max(i['N'] for i in iterations))
            tex.append(r'\expandafter\def\csname '+command+str(rho)+r'\endcsname{'+str(r['iterations'])+'}')
            coarse=next(r for r in iterations if r['contrast']==rho and r['method']==method and r['N']==min(i['N'] for i in iterations))
            tex.append(r'\expandafter\def\csname '+command+str(rho)+r'Coarse\endcsname{'+str(coarse['iterations'])+'}')
    tex.append(r'\newcommand{\ProfileN}{'+str(N)+'}')
    tex.append(r'\newcommand{\IterationN}{'+str(max(r['N'] for r in iterations))+'}')
    tex.append(r'\newcommand{\IterationCoarseN}{'+str(min(r['N'] for r in iterations))+'}')
    fdm_profile=read(f'layered_fdm_resistance_N{N}_profile')
    fem_profile=read(f'layered_fem_fitted_N{N}_profile')
    tex.append(r'\newcommand{\FDMSlice}{'+f"{fdm_profile['slice_y']:.3f}"+'}')
    tex.append(r'\newcommand{\FEMSlice}{'+f"{fem_profile['slice_y']:.3f}"+'}')
    (root/'slides'/'results_macros.tex').write_text('\n'.join(tex)+'\n',encoding='utf-8')
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
    fig,axs=plt.subplots(1,2,figsize=(9,4),sharex=True,sharey=True)
    iterations=read('iterations')
    for ax,method in zip(axs,('fdm_resistance','fem_fitted')):
        contrasts=sorted(set(r['contrast'] for r in iterations))
        cmap=plt.get_cmap('Blues' if method=='fdm_resistance' else 'Purples')
        shades=cmap(np.linspace(.45,.92,len(contrasts)))
        marker,line=('s','--') if method=='fdm_resistance' else ('x',':')
        for rho,shade in zip(contrasts,shades):
            rows=[r for r in iterations if r['method']==method and r['contrast']==rho]
            ax.loglog([r['N'] for r in rows],[r['iterations'] for r in rows],marker=marker,linestyle=line,
                      color=shade,linewidth=1.7,markersize=4.5,markerfacecolor='white',label=f'contrast {rho:g}')
        ns=sorted(set(r['N'] for r in iterations))
        ax.set_xticks(ns,labels=[str(n) for n in ns]);ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set(title=LABELS[method],xlabel='Intervals N');ax.legend(fontsize=9,frameon=False);ax.grid(alpha=.15,which='both')
    axs[0].set_ylabel('Unpreconditioned CG iterations')
    save(fig,'contrast_iterations')
    from .presentation_studies import figures
    figures(root)
