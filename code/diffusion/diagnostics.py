"""Physical norms, independent side fluxes, and algebraic-error correction."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import spsolve, norm as sparse_norm
from skfem import Basis, ElementTriP1


def fdm_l2(s, values, exact=True, interface_aware=True, order=4):
    """Integral L2 of piecewise bilinear reconstruction, split at interface.

    Resistance reconstruction uses accumulated resistance as the x coordinate
    on crossing cells. Ordinary bilinear reconstruction is also reported.
    For correction fields exact=False integrates the field itself.
    """
    z, w = leggauss(order); z=(z+1)/2; w=w/2
    total=0.
    y0=s.y[:-1][None,:]; dy=np.diff(s.y)[None,:]
    for i in range(s.N):
        lo, hi=s.x[i:i+2]
        cuts=[lo,hi]
        if s.model.name=='layered' and lo < .5 < hi:
            cuts=[lo,.5,hi]
        for low,high in zip(cuts[:-1],cuts[1:]):
            for zx,wx in zip(z,w):
                xx=low+(high-low)*zx
                t=(xx-lo)/(hi-lo)
                if interface_aware and s.model.name=='layered':
                    def resistance(end):
                        return (max(0.,min(end,.5)-lo) if lo<.5 else 0.) + max(0.,end-max(lo,.5))/s.model.contrast
                    t=resistance(xx)/resistance(hi)
                yy=y0+dy*z[:,None]
                lowvals=(1-t)*values[i,:-1]+t*values[i+1,:-1]
                highvals=(1-t)*values[i,1:]+t*values[i+1,1:]
                v=(1-z[:,None])*lowvals+z[:,None]*highvals
                e=v-s.model.exact(xx,yy) if exact else v
                total += (high-low)*wx*np.sum(w[:,None]*dy*e**2)
    return float(np.sqrt(total))


def physical_errors(s, values=None, exact=True):
    u=s.u if values is None else values
    if s.basis is None:
        return {'l2': fdm_l2(s,u,exact),
                'l2_bilinear': fdm_l2(s,u,exact,False)}
    basis=Basis(s.basis.mesh,ElementTriP1(),intorder=8)
    field=basis.interpolate(u); x=basis.global_coordinates()
    e=field-s.model.exact(*x) if exact else field
    grad=field.grad-s.model.gradient(*x) if exact else field.grad
    return {'l2':float(np.sqrt(np.sum(e**2*basis.dx))),
            'energy':float(np.sqrt(np.sum(s.model.a(*x)*np.sum(grad**2,axis=0)*basis.dx)))}


def metrics(s):
    u=s.u.ravel()[s.free]; r=s.b-s.A@u
    # One residual-correction solve estimates error to the exact discrete solve.
    # This is a numerical cross-check, not a certified roundoff bound.
    correction=np.zeros_like(s.u)
    correction.ravel()[s.free]=spsolve(s.A,r)
    alg=physical_errors(s,correction,exact=False)['l2']
    if s.basis is None:
        X,Y=np.meshgrid(s.x,s.y,indexing='ij'); ue=s.model.exact(X,Y)
        error=s.u-ue
        nodal_l2=np.linalg.norm(error[1:-1,1:-1])/s.N
    else:
        ue=s.model.exact(*s.basis.doflocs); error=s.u-ue; nodal_l2=None
    off=s.A.copy(); off.setdiag(0); off.eliminate_zeros()
    ans=dict(method=s.method,model=s.model.name,contrast=s.model.contrast,N=s.N,
             h=1/s.N,unknowns=len(s.b),nnz=s.A.nnz,
             assembly_s=s.assembly_s,solve_s=s.solve_s,
             nodal_linf=float(np.max(np.abs(error))),nodal_l2=nodal_l2,
             residual_inf=float(np.linalg.norm(r,np.inf)),
             residual_rel2=float(np.linalg.norm(r)/max(np.linalg.norm(s.b),np.finfo(float).tiny)),
             symmetry_rel=float(sparse_norm(s.A-s.A.T)/sparse_norm(s.A)),
             max_offdiagonal=float(max(0.,off.data.max(initial=0.))),
             algebraic_l2_correction=alg,minimum=float(s.u.min()),maximum=float(s.u.max()))
    ans.update(physical_errors(s))
    if s.model.name=='layered':
        ans['boundary_bound_violation']=float(max(0.,-s.u.min(),s.u.max()-1.))
        flux=flux_profile(s)
        qstar=-2*s.model.contrast/(1+s.model.contrast)
        ans['flux_relative_linf']=float(np.max(np.abs(flux['q']-qstar))/abs(qstar))
        ans['side_flux_jump_relative']=float(abs(flux['q_left']-flux['q_right'])/abs(qstar))
        ans['flux_slice_y']=flux['slice_y']
    return ans


def flux_profile(s):
    """x flux along actual central grid row / interior horizontal FEM slice."""
    if s.basis is None:
        j=s.N//2
        q=-s.ax[:,j]*np.diff(s.u[:,j])*s.N
        mid=(s.x[:-1]+s.x[1:])/2
        left=np.flatnonzero(s.x[1:]<=.5)[-1]
        right=np.flatnonzero(s.x[:-1]>=.5)[0]
        return dict(x=mid,q=q,temperature_x=s.x,temperature=s.u[:,j],
                    slice_y=float(s.y[j]),q_left=float(q[left]),q_right=float(q[right]))
    # Evaluate gradients in elements intersected by a non-grid-line slice.
    x=np.sort(np.r_[np.linspace(.001,.499,160),np.linspace(.501,.999,160)])
    yy=np.full_like(x,.5+.173/s.N)
    elements=s.basis.mesh.element_finder()(x,yy)
    field=s.basis.interpolate(s.u)
    q=-s.model.a(x,yy)*field.grad[0,elements,0]
    xx=np.r_[0,x,1.]; ys=np.full_like(xx,yy[0])
    temp=s.basis.interpolator(s.u)(np.array([xx,ys]))
    return dict(x=x,q=q,temperature_x=xx,temperature=temp,slice_y=float(yy[0]),
                q_left=float(q[np.flatnonzero(x<.5)[-1]]),
                q_right=float(q[np.flatnonzero(x>.5)[0]]))
