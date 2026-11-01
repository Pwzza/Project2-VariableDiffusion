"""Independent checks of assembly, nonzero boundaries, flux and convergence.

Run: .venv/Scripts/python.exe -m unittest discover -s tests -v
"""
from pathlib import Path
import sys
import unittest
import numpy as np
from scipy.sparse.linalg import eigsh
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from diffusion.models import Model
from diffusion.solvers import assemble_fdm, assemble_fem, face_coefficient, solve
from diffusion.diagnostics import metrics, physical_errors, flux_profile
from diffusion.experiments import independent_forcing_check


class NumericalChecks(unittest.TestCase):
    def test_constant_coefficient_hand_matrix_both_methods(self):
        expected=np.array([[4,-1,-1,0],[-1,4,0,-1],[-1,0,4,-1],[0,-1,-1,4]])
        fd=assemble_fdm(Model('layered',1),3,'resistance')
        # FEM with rho=1 still inserts the interface, so use the ordinary
        # smooth mesh but a test-only constant coefficient model.
        class Constant(Model):
            def a(self,x,y):
                return np.ones(np.broadcast_shapes(np.shape(x),np.shape(y)))
        fe=assemble_fem(Constant(),3)
        np.testing.assert_allclose(fd.A.toarray()/9,expected,atol=1e-13)
        np.testing.assert_allclose(fe.A.toarray(),expected,atol=1e-13)

    def test_forcing_independent_flux_differentiation(self):
        self.assertLess(independent_forcing_check()['maximum_absolute_defect'],1e-7)

    def test_hand_fdm_row(self):
        s=assemble_fdm(Model(),4)
        # Point (.25,.25): shared x-face coefficients 1.125,1.375,
        # y-face coefficients 1.25,1.25; h^-2=16.
        self.assertAlmostEqual(s.A[0,0],80.)
        self.assertAlmostEqual(s.A[0,3],-22.)
        self.assertAlmostEqual(s.A[0,1],-20.)
        self.assertAlmostEqual(s.b[0],float(Model().forcing(.25,.25)))

    def test_general_interface_resistance(self):
        m=Model('layered',10)
        a=face_coefficient(m,.4,.7,0.,0.,'resistance')
        self.assertAlmostEqual(float(a),.3/(.1+.2/10))

    def test_boundary_elimination_and_layered_exactness(self):
        for rho in (.1,1.,10.,100.):
            for constructor in (lambda:assemble_fdm(Model('layered',rho),15,'resistance'),
                                lambda:assemble_fem(Model('layered',rho),15)):
                s=solve(constructor());r=metrics(s)
                self.assertLess(r['nodal_linf'],1e-11)
                self.assertLess(r['flux_relative_linf'],1e-10)
                self.assertLess(r['boundary_bound_violation'],1e-12)
                self.assertLess(r['symmetry_rel'],1e-13)
                self.assertGreater(eigsh(s.A,k=1,which='SA',return_eigenvectors=False)[0],0)
                self.assertLess(r['max_offdiagonal'],1e-12)

    def test_arithmetic_is_symmetric_but_physically_wrong(self):
        s=solve(assemble_fdm(Model('layered'),15,'arithmetic'));r=metrics(s)
        self.assertLess(r['symmetry_rel'],1e-13)
        self.assertGreater(r['nodal_linf'],1e-3)
        self.assertGreater(r['flux_relative_linf'],1e-2)
        self.assertLess(r['boundary_bound_violation'],1e-12)

    def test_reconstruction_distinguishes_exact_nodes_from_kink(self):
        s=solve(assemble_fdm(Model('layered'),15,'resistance'))
        e=physical_errors(s)
        self.assertLess(e['l2'],1e-12)
        self.assertGreater(e['l2_bilinear'],1e-4)

    def test_smooth_convergence_and_algebraic_separation(self):
        for assemble in (assemble_fdm,assemble_fem):
            results=[metrics(solve(assemble(Model(),N))) for N in (8,16,32)]
            p=np.log(results[-2]['l2']/results[-1]['l2'])/np.log(2)
            self.assertTrue(1.8<p<2.2,p)
            for r in results:
                self.assertLess(r['algebraic_l2_correction'],1e-3*r['l2'])

    def test_fitted_mesh_does_not_cross_interface(self):
        s=assemble_fem(Model('layered'),15)
        vertices=s.basis.mesh.p[0,s.basis.mesh.t]
        self.assertFalse(np.any((vertices.min(axis=0)<.5)&(vertices.max(axis=0)>.5)))


if __name__=='__main__':
    unittest.main()
