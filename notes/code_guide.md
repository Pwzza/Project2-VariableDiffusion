# Numerical implementation and presentation evidence

Current implementation: 2026-11-01. This guide supersedes the implementation-status and experiment-numbering statements in the original October design notes. Preserve those notes as pre-computation history, not as the current code specification.

CG study update (2026-11-08): full mode now uses N=15,31,63,127 and contrasts 1,3,10,30,100 for both methods (40 cases). Run `python run_all.py --iterations-only` to update this study and its plots without rerunning the other experiments. The recomputed final residual must be at most 1e-10; CG requests 9e-11 to allow a small margin. Every row saves termination status, unknown count and nodal maximum error. The two-panel logarithmic plot keeps resistance FDM blue and fitted FEM purple, distinguishing contrasts with markers/line styles. The slide presents the observations and leaves “Why?” open. Original iteration-study settings below describe the earlier implementation.

Presentation update: a nine-frame Madrid Beamer draft now exists in slides/main.tex, with a compiled preview in slides/main.pdf. Timing is left to the presenter. Plot regeneration also exports the exact-flux short-vector map and selected slide figures, and generates slides/results_macros.tex directly from saved results. The numerical solves are unchanged.

## Story and scope

The presenter manages timing. The nine-slide draft emphasizes interface physics and how FDM/FEM preserve it. There is no standalone Experiment 0: preliminary analytic visuals introduce both models and the two methods. Experiment 1 compares both methods on smooth diffusion. Experiment 2 compares both methods on the layered model, including the deliberately incorrect arithmetic FDM closure. Error diagnostics sit inside each experiment; optional cost analysis follows the method comparisons.

## Modules

| File in code/diffusion | Responsibility |
|---|---|
| models.py | Immutable smooth/layered benchmark, analytic gradient/flux, full-boundary g, variable contrast |
| solvers.py | COO-to-CSR shared-face FDM; scikit-fem 2D P1 stiffness/load; strong boundary elimination; direct reference solve |
| diagnostics.py | Physical and nodal norms, residual correction, sign/symmetry checks, side-flux reconstruction |
| experiments.py | Four meshes, quadrature check, eigenvalue corroboration, CG contrast study, saved data, acceptance gates |
| plots.py | Headless mathematical introduction/method diagrams and figures from saved numerical data |

Provenance: the design follows Materials/03_Scaffolding/poisson_fd_skeleton.py, fem2d_scikitfem_skeleton.py, convergence_plot_template.py and iterative_solvers_skeleton.py. These are course references, not instructions overriding the approved plan. The original course setup checks remain unchanged in code/course_checks. Project solver modules are new implementations, not completed TODO edits to the supplied files.

## Discretization conventions

Domain is the unit square. N always means number of macro intervals, not number of interior unknowns. h=1/N is the common macro-grid spacing; FEM triangle diameter on square cells is sqrt(2)h. Actual unknown counts are saved.

FDM arrays use U[i,j]=u(x_i,y_j), with j changing fastest. The reduced interior index is (i-1)(N-1)+(j-1). A single positive face coefficient is shared between adjacent rows. Boundary contributions move into b; Dirichlet rows are not overwritten with identity. Smooth coefficients use arithmetic midpoint values, exact for a=1+x. Layered resistance coefficients integrate the subedge resistances; arithmetic is a deliberate comparison.

FEM uses right triangles from MeshTri.init_tensor and P1 basis functions. The layered mesh inserts x=0.5 into the same macro grid, creating a fitted interface and extra degrees of freedom. No element straddles the jump. Quadrature order 6 assembles stiffness/load; order 8 integrates errors, and a separate finest-mesh assembly at order 8 checks quadrature sensitivity. All boundary degrees of freedom receive exact g, and are eliminated from the linear system.

## Experiments

1. Smooth: N=16,32,64,128; arithmetic FDM and P1 FEM. Fit the last three levels with saved window and maximum log deviation. Report continuous L2, nodal maximum, FDM scaled nodal L2, FEM energy, residuals and residual-correction estimates. No universal FEM maximum-norm rate is assumed.
2. Layered: N=31,63,127,255; arithmetic FDM, resistance FDM, fitted FEM. Exact forcing is zero and exact g is applied on all four sides. Save temperature/flux profiles, error norms, independent side flux mismatch and boundary bounds. Do not fit roundoff-level errors.
3. Computational difficulty within the layered study: contrasts 1,10,100 and N=15,31,63, both methods. Unpreconditioned CG, zero initial guess, relative Euclidean residual target 1e-10, absolute target zero. Fail on nonconvergence or inadequate final residual. Counts are not directly comparable to an unspecified preconditioned solver.

Quick mode uses four smaller meshes and a smaller iteration study; metadata.json distinguishes it from the final configuration. Quick mode overwrites the same filenames; stale case files from other configurations can remain. Use the N values listed in metadata and tables, not a directory glob, to identify current results. Plot generation follows those tables.

## Norms and interpretation

- Continuous L2 is an integral over the domain. Smooth FDM uses bilinear reconstruction integrated with tensor Gauss quadrature. FEM uses integrated P1 error.
- Layered FDM's primary L2 uses resistance-coordinate reconstruction on crossing cells, split at the interface for integration. Both arithmetic and resistance solutions use this explicitly stated reconstruction. Ordinary bilinear reconstruction error is separately saved as l2_bilinear. Exact nodal values alone do not make a bilinear interpolant exact across a kink.
- nodal_linf is the maximum error at degrees of freedom; it is not a quadrature-based continuous supremum. FDM nodal_l2=h times the interior Euclidean error norm. FEM energy integrates a times squared gradient error.
- residual_rel2=||b-AU||2/||b||2 on the reduced system. Matrix scaling differs between FDM and FEM, so raw residuals are not physical accuracy comparisons.
- algebraic_l2_correction is the physical norm of one residual-correction solve A delta=b-AU. It is a numerical estimate/cross-check, not a rigorous floating-point error bound. Smooth acceptance requires it below 0.001 times the discretization error. Exactly representable layered cases instead use absolute nodal and relative flux thresholds.
- Flux errors are maxima along the sampled horizontal slice, not whole-domain flux norms. Actual slice y values are saved and labelled. Odd FDM N has no y=0.5 row. FEM uses a slightly offset interior slice to avoid ambiguous triangle-edge evaluation.
- FDM flux uses shared face differences. Independent side mismatch compares adjacent same-material edges on either side of the crossing edge; it does not reuse the single interface face twice. FEM uses element gradients on each side. An arithmetic solution can have small side mismatch but incorrect flux magnitude, so both diagnostics are necessary.
- Boundary bounds and a nonpositive off-diagonal sign pattern support the discrete maximum principle on these meshes; they do not prove interface accuracy. SPD follows from the positive discrete/weak energy after elimination. Smallest eigenvalues on a representative mesh are numerical corroboration only.

## Output map

data/smooth.csv and layered.csv contain all mesh metrics; JSON preserves explicit null/missing norms. fits.json contains selected fit windows. checks.json contains forcing, quadrature, representative eigenvalue checks and acceptance gates. iterations.csv contains counts, attained residuals and timings. metadata.json records package versions, actual mesh lists and run mode. Per-case NPZ files contain nodal fields and mesh geometry, and layered profile JSON contains sample coordinates and fluxes.

Figures export PDF and PNG: intro_models, intro_interface, methods_geometry, fem_basis, interface_resistance, matrix_structure, smooth_convergence, smooth_error_maps, layered_profiles, layered_errors, contrast_iterations. Analytic visuals are not empirical evidence. Result plots read saved tables/fields. Exact layered curves can overlap; layered_errors and the temperature-error panel expose differences.

Single-run timings are diagnostic and machine-dependent. Controlled cost-versus-accuracy benchmarks, AMG comparisons, unresolved FEM, and the final slide assembly are not implemented. No speed ranking should be claimed from the present timings. Theoretical assembly/storage is O(h^-2); solver complexity needs its own ordering/preconditioning assumptions.

## Reproduction and validation

From NumericalPDE, run `.venv/Scripts/python.exe -m unittest discover -s tests -v`, then `.venv/Scripts/python.exe run_all.py`. Paths resolve relative to run_all.py, so it also works from another current directory. No notebook state or interactive plotting is required.

Tests cover a hand-derived 4-by-4 constant-coefficient matrix for both methods, a variable-coefficient boundary row, arbitrary-distance resistance, independent forcing differentiation, full-boundary reproduction across several contrasts, matrix positivity/symmetry/sign, fitted mesh alignment, arithmetic closure failure, reconstruction at a kink, smooth convergence and algebraic separation. These checks are distinct from the supplied environment tests.

Human responsibilities remain: independent review of formulas and solver code, reproducing outputs, course audit and individual-contribution records. Automated checks are not a substitute for team sign-off.
