# Experiment design and theoretical predictions

**Proposed title: Reliable heat flow through heterogeneous materials: accuracy, interface physics, and computational cost.**

This document plans experiments; no predicted slope, timing, error, or solver advantage below is a measured project result. Complete the prediction/output/reconciliation record before accepting each future experiment.

## 1. The story connecting the experiments

Project 1 separated exact-reference convergence, structural stress testing, numerical-reference reliability, and cost. Retain that logic here, but replace stochastic sampling error with algebraic solver error, and positivity of stochastic paths with conservation and the discrete maximum principle.

Our practical question is: *If a composite has a low-conductivity layer, can we trust both the temperature and the heat flow predicted by our solver?*

The chain of questions is:

1. **Can we trust the assembly?** E0 tests signs, indexing, nonzero boundary data and simple known limits.
2. **Does refinement buy the accuracy theory promises?** E1 establishes smooth second-order accuracy against exact data.
3. **Does smooth success survive a material jump?** E2 exposes a wrong-resistance closure that still looks stable and symmetric.
4. **Are we seeing PDE error or an unfinished linear solve?** E3 measures the algebraic error floor.
5. **Which solver is worthwhile at the accuracy we actually need?** E4 is the substantial depth comparison, including refinement, contrast and cost.
6. **Can the good interface rule handle geometry beyond a midpoint?** Optional E5 generalizes to an off-center interface using exact resistance.

The central claim to test is conditional: a conservative face scheme with correct resistance should deliver smooth accuracy, correct interface transport, and a defensible cost-to-accuracy tradeoff. None of these is certified by a smooth plot or a small residual alone.

## 2. Fixed definitions and shared measurement protocol

Use the conventions in `problem_definition.md`. For nodal error `e=U-u*`, always report

\[
E_{2,h}=\left(h^2\sum_{i,j=1}^{N-1}e_{ij}^2\right)^{1/2}=h\|e\|_2,
\qquad E_\infty=\max_{i,j}|e_{ij}|.
\]

These are a discrete interior L2 norm and a nodal maximum norm, not an unqualified continuous L2 integral. Optional discrete energy error is `E_A=sqrt(h^2 e^T A e)`, defined on the nodal comparison error with zero boundary error; do not equate it automatically to a continuous gradient error.

Report the truncation residual `tau=b-A*u*_nodes` separately from algebraic residual `r=b-A*U`. The correct identity with this definition is `A*(U-u*_nodes)=tau` for an exact discrete solve. Record `h*||tau||2` and `||tau||inf`; a residual is converted to error by `A^{-1}`, not merely by changing its norm weights.

Solver diagnostics on every row of every refinement table:

- `r_abs_inf=||r||inf` and `r_rel_2=||r||2/||b||2` (if b=0, report absolute residual; do not divide by zero).
- Backward error `eta_inf=||r||inf/(||A||inf*||U||inf+||b||inf)`.
- Solver name, requested tolerance or `N/A (direct)`, termination flag, iterations and iteration cap.
- Relative asymmetry `||A-A.T||F/||A||F`, DOFs, nnz and wall-time components.

Use two-grid rates `p=log(E_old/E_new)/log(h_old/h_new)` and a least-squares fit of `log(E)` against `log(h)` on at least three demonstrably asymptotic levels. Show all four levels, pairwise rates, selected fit levels, RMS log-fit residual and sensitivity of the fitted slope to dropping the coarsest point. No artificial slope fit to round-off-level data. A practical review flag is a smooth finest-three fit outside 1.8-2.2; investigate rather than altering data to meet it.

Benchmark timing: one warm-up; five measured repetitions (three for expensive cases); report median and interquartile range. Use identical matrix/RHS, machine, initial guess (zero), stopping rule and thread configuration. Record Python/package versions, CPU/OS, mesh, contrast, command, revision, and run parameters. Separate assembly, solver setup/factorization/preconditioner setup, solve, diagnostics and plotting. Include setup in a single-RHS total; show reused setup separately only for an explicitly stated multiple-RHS scenario. Extra residual logging is diagnostic overhead: record its sampling frequency and treat methods consistently.

## E0. Establish trust before interesting physics

**Motivation:** A transposed image, missing Dirichlet term or sign error can survive a homogeneous, symmetric test.

**Inputs and checks:**

1. Supplied complete Poisson setup checker and analytic-bank self-check. These prove environment readiness and input consistency, not our future solver.
2. Our assembly at `N=3,a=1`: match all entries of the hand matrix in the definition sheet; verify map/inverse map, CSR, sign pattern and expected `nnz=5*n^2-4*n` for `n=N-1`.
3. `a=1,u*=x^2-y^2,f=0,g=u*`: exact to round-off because central differences reproduce quadratics. Exercises nonzero data on all sides. Do not fit an order.
4. Tutorial P5, `a=1,u*=exp(x)sin(pi*y)`, `f=(pi^2-1)u*`. At `N=16` predict the center near `1.652` from the bank; exact PDE value is `exp(1/2)`, not that discrete number. Repeat `N=32` to exercise boundary rows and reduce error. This asymmetric solution also catches x/y plotting swaps.
5. Constant-coefficient scaling: for fixed `f,g=0`, replacing `a` by `c*a` scales `U` by `1/c`. For `f=0` and fixed boundary data, common scaling of all conductivities leaves U unchanged and scales flux by c.

**Expected evidence:** small matrix and two boundary-adjacent rows; bank output; values/errors/scaling defects in `e0_checks.csv`; one spy plot and one deliberately asymmetric heatmap.

**Gate:** round-off-scale defects assessed with matrix/solution scale and conditioning; no fixed absolute residual target imposed across arbitrary mesh/contrast. Stop if the hand matrix, signs, boundary shift or indexing is wrong.

**Bridge:** Correct simple cases establish implementation basics; they cannot establish an asymptotic accuracy rate for heterogeneous diffusion.

## E1. Verify smooth variable-coefficient accuracy

**Motivation:** A code for `-a*Delta(u)` solves the wrong PDE yet can pass a constant-coefficient test. Force the coefficient-gradient term to matter.

Use the prescribed smooth benchmark:

\[
a=1+x,\quad u_*=\sin(\pi x)\sin(\pi y),\quad g=0,
\]
\[
f=-a\Delta u_*-\nabla a\cdot\nabla u_*
=2\pi^2(1+x)\sin(\pi x)\sin(\pi y)
-\pi\cos(\pi x)\sin(\pi y).
\]

**Independent input verification:** manually differentiate `u_x,u_xx,u_yy`; compare this expression to the course bank P4 and its self-check. Separately approximate `-div(a grad u*)` by centered face fluxes at several off-grid interior points while reducing the finite-difference step, looking for second-order decrease before round-off. This checker must not call project assembly or simply reuse its forcing function.

**Theory:** In one dimension the arithmetic face average satisfies

\[
\frac{a(x)+a(x+h)}2=a(x+h/2)+\frac{h^2}{8}a''(x+h/2)+O(h^4).
\]

The centered edge derivative is `u'(x+h/2)+h^2*u'''(x+h/2)/24+O(h^4)`. With smooth a and u, taking the centered difference of the two edge fluxes yields second-order local consistency. For `a=1+x` the arithmetic coefficient is exactly the midpoint value. Uniform ellipticity, the discrete comparison principle and a smooth barrier provide stability; under the smooth square-domain data here this gives expected `Einf=O(h^2)` and hence `E2h=O(h^2)` (the latter is bounded by the maximum error on this unit-area grid). State the regularity assumptions: this argument is not a theorem for discontinuous coefficients or arbitrary domains.

For the one-dimensional part with linear a, the leading discrete-minus-continuous operator defect is

\[
(L_h-L)u=-h^2\left(\frac{a\,u^{(4)}}{12}+\frac{a'\,u^{(3)}}6\right)+O(h^4).
\]

Sum the x and y contributions. The sign reverses for `tau=L(u)-L_h(u)`; it is its magnitude and order that we fit.

**Protocol:** `N=16,32,64,128` (exact `h,h/2,h/4,h/8`), then `256` only to resolve a rate concern or extend cost measurement. Primary baseline is `spsolve`, with direct residual/forward-error evidence described in E3. Plot exact, numerical and signed-error surfaces; log-log errors and truncation residuals with reference `h^2` slopes.

**Qualitative predictions:** u is nonnegative and symmetric about `y=1/2`; u* is also symmetric about `x=1/2`, even though a and f individually are not. The discrete problem enforces y symmetry, while the x symmetry is only exact-solution structure and may have a discretization-scale defect. Do not claim that f is everywhere nonnegative: it changes sign near the left boundary. The zero-boundary maximum-principle bound for harmonic problems does not apply to this nonzero-source case.

**Falsification:** in a separate rejected variant omit `-a_x*u_x` from f while keeping the intended exact solution. Predict a nonvanishing solution-error floor against u*, despite a tiny algebraic residual. Restore the proper f for all accepted results. This demonstrates that solver convergence cannot validate PDE transcription.

**Optional stronger nonzero-boundary MMS:** `a=1+x`, `u*=exp(x)sin(pi*y)`,
`f=(((1+x)*(pi^2-1))-1)*u*`, `g=u*`. Run at least two levels; use it to test variable coefficients and boundary elimination together. It is stronger than merely repeating a zero-boundary case.

**Bridge:** Smooth second order validates the method where derivatives exist. A composite introduces a gradient jump, so we next test the physical quantity the model actually conserves.

## E2. Interface stress test: stability can coexist with wrong physics

**Motivation:** Arithmetic and harmonic face means are both positive, symmetric, monotone choices. Only one represents resistances in series across a jump.

Use `aL=1,aR=10`, interface `xi=1/2`, `f=0`. Derive the exact solution from constant `a*u_x=C`, continuity and the end temperatures 0 and 1:

\[
R=\frac{\xi}{a_L}+\frac{1-\xi}{a_R}=\frac{11}{20},
\quad C=1/R=20/11,\quad q_x=-C=-20/11,
\]
\[
u_*(x)=\begin{cases}(20/11)x,&x\le1/2,\\(2x+9)/11,&x\ge1/2.\end{cases}
\]

Set `g=u*(x)` on **every** boundary edge. Thus this x-only exact solution is also a valid 2D Dirichlet benchmark. Continuity gives interface temperature `10/11`, and left/right slopes differ by a factor of 10.

Across an edge containing portions dL and dR, integrating `u_x=C/a` yields

\[
U_R-U_L=C(d_L/a_L+d_R/a_R),\qquad
q_{LR}=-a_{\rm face}(U_R-U_L)/h.
\]

This derives the resistance coefficient. For a bisected edge `aH=20/11`, whereas `aA=11/2`. Inserting the *exact nodal profile* into the arithmetic crossing flux gives

\[
\frac{q_A[u_*]}{q_*}=\frac{a_A}{a_H}
=\frac{(a_L+a_R)^2}{4a_La_R}=\frac{121}{40}=3.025.
\]

That ratio is an analytic **consistency test**, not the flux ratio of the re-solved arithmetic system. Re-solving changes U. The arithmetic numerical solution can gain y-dependence because the top/bottom Dirichlet values still impose the correct exact profile. Consequently a single 1D total-resistance formula does not describe that entire wrong-closure 2D solution.

**Protocol:** `N=31,63,127,255` intervals as prescribed: each crossing edge is bisected and no interface node exists. Use actual `h=1/N` in rates; this sequence is not exact halving. Use the same grid, f, g and solver for the two face rules. Plot a clearly labelled slice on a grid row close to `y=1/2` (odd N has no center row), its neighboring face fluxes, a 2D error map, and temperature/flux errors versus h.

**Predictions:**

- Harmonic closure reproduces exact nodal u and constant x flux to linear-solve/round-off accuracy. No meaningful convergence slope can be extracted from this error floor.
- Arithmetic closure fails exact reproduction and produces an interface error even though SPD and boundary bounds remain valid. An approximately first-order solution-error trend is a hypothesis to measure for this 2D setup, not an unconditional theorem. The interface flux can behave differently from the global solution norm; measure it separately.
- Both solutions should remain within `[0,1]`, up to algebraic and floating-point error. This confirms the maximum-principle structure and does not decide which temperature field is accurate.

**Flux diagnostics that cannot be replaced by a zero jump:**

1. Numerical face flux `Qx=-a_face*(U[i+1,j]-U[i,j])/h`, with exact signed target `-20/11`; report relative crossing-face flux error across all interior y rows and on the central slice.
2. Algebraic conservation: neighboring rows use opposite signed copies of the **same face flux**. Their internal cancellation is exact by construction for both schemes. Report it as an assembly check, never as independent validation of the averaging rule.
3. Independently reconstructed material-side fluxes: use only same-material nodes on each side to obtain second-order one-sided derivatives at the interface (fit a quadratic to three distinct columns and differentiate it at `xi`); form `qL=-aL*uxL`, `qR=-aR*uxR` using the same +x orientation. Record `max_j|qL-qR|/|q*|` and each side's error relative to q*. Avoid fitting across the kink. For the exact piecewise-linear profile this reconstruction is exact up to round-off.
4. Report the jump between crossing-face flux and the nearest same-material x-face fluxes, and the full cell divergence of x and y fluxes. Those neighboring x fluxes need not match for a y-dependent numerical solution; distinguish that diagnostic from a genuine interface transmission jump.
5. Boundary-integrated outward flux versus integrated source, with a declared boundary-segment quadrature and corner convention. A finite-difference nodal control-volume conservation test must use its own represented control-volume area, not casually equate an interior-node sum with the full domain integral.

**Bridge:** The exact interface case reveals physical mistakes, but its harmonic solution has no useful discretization-error scale. Use the smooth benchmark to quantify how tightly the solver must be run before evaluating efficiency.

## E3. Separate algebraic error from PDE error

**Motivation:** Refinement increases conditioning sensitivity. A nominally small residual can still obscure a second-order convergence table, particularly with nonzero boundary contributions making b large.

Let `U_D` be the direct discrete solve, `U_k` an iterate and `u*_h` the exact PDE nodal values. Decompose

\[
U_k-u_*^h=(U_D-u_*^h)+(U_k-U_D).
\]

The first term is discretization error (up to direct round-off); the second is algebraic error relative to the discrete reference. They can cancel, so comparing total errors alone is insufficient.

For SPD A, `||U_k-U_D||2 <= ||r_k||2/lambda_min(A)` when U_D denotes the exact system solution. For approximate direct reference include its own residual bound. Since each edge coefficient is at least amin,

\[
\lambda_{\min}(A)\ge a_{\min}\,8h^{-2}\sin^2(\pi/(2N)).
\]

This follows from the energy comparison to the constant-coefficient five-point Dirichlet matrix. It supplies a conservative **absolute** bound in our measured norm:

\[
E_{\mathrm{alg},2h}\le h\|r\|_2/lambda_{\mathrm{lower}}.
\]

Use it for both direct and iterative residuals. On a moderate mesh corroborate eigenvalue estimates using `eigsh`, report convergence of the eigensolver, and do not present a numerical smallest-eigenvalue estimate as a certified lower bound. The familiar relative-error estimate `relative forward error <= kappa2(A)*relative residual` also explains why the residual alone is insufficient.

**Protocol:** on E1's finest fitted mesh (`N=128` initially), run residual-stopped Gauss-Seidel at `1e-4,1e-6,1e-8,1e-10,1e-12`. Log success/cap, true residual, E2h, Einf and an area-average QoI `h^2*sum(U)`; for the exact smooth solution the continuous integral is `4/pi^2`, but this discrete sum has quadrature error too. Compare against U_D in the same norms. Repeat a compact sweep on each remaining refinement level; tighten on contrast cases as needed.

**Acceptance:** algebraic E2h error or its residual-to-error bound must be at most `1e-3` of direct E1 discretization E2h error on every fitted level, especially the finest. For maximum-norm claims, independently require the analogous max-error comparison to U_D; a conservative bound is `||e||inf <= ||e||2`. Establish the reference via its residual bound and an independent tight solve; if using only successive solve-pairs, tighten twice and require the second difference to be <=10% of the first before treating the tighter solve as adequate. Do not demand an unattainable `1e-16` tolerance: stop with an explicit failed gate if round-off or iteration limits prevent the required evidence.

For the layered exact-reproduction case, discretization error is theoretically zero; a ratio to that error is meaningless. Instead require `Einf <= 1e-9` and relative flux error `<=1e-8` initially, adjusted only with documented conditioning/round-off evidence. These are proposed numerical verification targets, not pre-observed achievements.

**Output:** tolerance/error/residual/QoI table and error-floor plot. Select the loosest tested tolerance that meets the gates. Direct solver has no requested tolerance: record N/A, residual, backward error and forward-error bound.

**Bridge:** Once solver error is controlled, timing comparisons answer a real question: the work required for an equally trustworthy numerical answer.

## E4. Selected depth: direct versus iterative at matched accuracy

This is the explicit **Depth direction 4.2** comparison. Keep it distinct from Direction 3's mandatory arithmetic/harmonic comparison. Use `spsolve` versus Gauss-Seidel as the two required variants; CG and AMG-CG are optional tutorial-inspired third variants.

**Motivation:** Sparse direct solves are robust for moderate meshes, but refinement and material contrast can make computational cost dominate. We want a practical choice based on verified accuracy, not fewer iterations in isolation.

**Grid study:** smooth E1 data, `N=16,32,64,128` and optionally 256 if feasible. One direct solution per system outside timing repeats; timed solves may repeat on the same immutable matrix/RHS. Apply E3's residual rule and error gate. Report unknowns, nnz, iterations, assembly/setup/solve/total time, error and memory indicators. Sparse matrix bytes and factor nnz are legitimate storage measures; label them separately from measured process peak memory.

**Contrast study:** fix `N=63` then confirm selected cases at `N=127`. Run the layered exact family with `rho=aR/aL=1,10,100,1000`, `aL=1` and

\[
C(\rho)=\frac{2\rho}{1+\rho},\qquad
u_*(x)=\begin{cases}C x,&x\le1/2,\\C[1/2+(x-1/2)/\rho],&x\ge1/2,\end{cases}
\quad q_*=-C.
\]

Generate the corresponding g afresh on all four sides. Use harmonic faces for the valid comparison. Because the exact nodal solution is available, report actual nodal algebraic error as well as residual; normalized residuals alone can be misleading as boundary terms/contrast scale. Use a common nominal tolerance initially, then tighten cases to meet the same E3 layered temperature/flux targets, reporting every chosen tolerance.

**Theory:** energy comparison gives

\[
a_{\min}A_0\preceq A\preceq a_{\max}A_0,\qquad
\kappa_2(A)\le\frac{a_{\max}}{a_{\min}}\kappa_2(A_0)
=\rho\,O(h^{-2}).
\]

This is a bound, not an exact contrast law or a guarantee of monotonically increasing observed iterations. Uniform scaling of all a changes eigenvalues but not kappa. Estimate extremal eigenvalues only on affordable meshes, separately from solver timing.

For constant-coefficient regular grids, Jacobi has spectral radius `cos(pi*h)=1-O(h^2)`; lexicographic GS is faster but still needs order `h^-2*log(1/tol)` sweeps in the standard setting. A sweep costs O(m) for `m=(N-1)^2`, so stationary methods become expensive under refinement. Heterogeneous cases need measured rates; do not reuse the Poisson-optimal SOR omega without evidence. Generic unpreconditioned CG has worst-case iteration bound proportional to `sqrt(kappa)*log(1/tol)`, not an exact iteration prediction. Near mesh-independent AMG-CG convergence is a hypothesis conditional on good coarsening and coefficient handling.

Sparse direct storage/time depend on ordering and fill. The common 2D nested-dissection estimates O(m log m) storage and O(m^(3/2)) factorization work are context, not measured guarantees for SciPy's actual ordering. Report ordering/backend and factor fill if available; never claim a universal crossover from asymptotic estimates.

**Matched-accuracy plots:** (i) time vs E2h on E1 across meshes, with target `E2h<=1e-3` and, if reached, `1e-4`; (ii) iterations vs h; (iii) iterations/condition estimates vs rho; (iv) residual histories on representative cases. Compare cost-to-target only among runs that actually pass tolerance/error gates. Do not interpolate a crossing unsupported by observed meshes.

Use a declared cap (initially 50,000 GS sweeps, or 120 s per solve, with checks every 10 sweeps); record cap hits as failures, not converged timing points. Count sweeps separately from residual checks. Implement GS via triangular solves or justified compiled loops; a slow Python node loop would confound algorithm and language costs. Run cost measurements serially on this machine.

**Expected conclusion form:** direct is likely practical for these moderate systems; GS likely worsens with refinement. An efficiency crossover or AMG benefit must be measured. If no crossover occurs, that is the result. Do not promise that iterative wins.

## E5. Optional geometry extension: exact resistance away from midpoints

Only after E0-E4 pass. Set interface `xi=0.37`, `aL=1,aR=10`, fixed uniform grid families avoiding interface nodes. Derive

\[
R=\xi/a_L+(1-\xi)/a_R,\quad C=1/R,\quad
u_*(x)=C\begin{cases}x/a_L,&x\le\xi,\\\xi/a_L+(x-\xi)/a_R,&x\ge\xi.\end{cases}
\]

Use the exact subedge distances for every crossing face. Compare distance-weighted resistance to an incorrectly assumed 50/50 harmonic mean. The former should again reproduce the nodal profile up to solve error; the latter misrepresents the interface's position. This is a low-cost, physically motivated extension of E2, with no need for another PDE topic. Report actual crossing fractions on every mesh.

## 3. Other depth directions considered and deferred

- **Two boundary closures:** would require adding Neumann/Robin conditions, changing the present all-Dirichlet study. Defer to preserve scope.
- **FDM vs 2D P1 FEM:** valuable independent method check using scikit-fem, and a good substitute for E4 if team interests change. Match the PDE, mesh scale, DoFs, algebraic tolerance and **continuous L2 error via common quadrature of declared reconstructions**, with nodal Linf on common nodes. On smooth regular data expect FEM L2 order 2 and H1/energy order 1. The FDM discrete nodal L2 error and FEM integrated L2 error are not automatically the same metric. For the layered FEM case use an interface-fitted mesh; odd-interval FDM and fitted FEM are different mesh families and must be disclosed. One-dimensional FEM is not an eligible 2D method comparison.
- **Smooth vs singular data / L-shaped domain:** would introduce domain/regularity effects separate from the chosen interface issue. The pack's FEM energy exponent about 2/3 and L2 exponent about 4/3 concern that singular-domain FEM problem; do not transfer them to our square-interface FDM data. Defer along with adaptive refinement, pure-Neumann nullspaces and nonlinear MEMS.
- **CG/AMG:** good extension inspired by Tutorial 6, especially as the required packages are installed. Include preconditioner setup and verify symmetric positive-definite preconditioning before CG. Do not silently substitute CG for the selected GS depth study.

## 4. Definition of a successful finished project

The final report's conclusions must be supported by E0-E4 actual outputs, at least two independent checks (analytic references plus a separate solver/hand derivation), a named four-mesh convergence study, nonzero-boundary evidence, an interface solution-and-flux comparison, matrix-structure evidence, and algebraic-error separation. Save predictions before runs and reconcile failed expectations openly. E5 is optional. The required report is at least 8 LaTeX pages, the presentation is 10 minutes plus 3 minutes Q&A, and the code archive must reproduce figures through `run_all.py`.
