# What to build next

Status: implementation specification. Only the environment driver and supplied setup checks currently exist. Do not confuse the module names below with delivered numerical solvers.

## Build order and interfaces

| Stage | Proposed module | Responsibility and review gate |
|---|---|---|
| 1 | `code/problems.py` | Immutable problem definitions: coefficient, forcing, exact u, boundary values, exact flux where known, interface position. Record N as intervals. Keep independently derived oracle formulas separate. |
| 2 | `code/grid.py` | Full/interior axes, ij grids, index/inverse map, boundary reconstruction. Check the N=3 matrix ordering and asymmetric plots. |
| 3 | `code/faces.py` | Construct full horizontal arrays of shape `(N,N+1)` and vertical arrays `(N+1,N)`; each edge stored once. Smooth arithmetic mean and jump arithmetic/resistance variants. Explicit distance-weighted formula for E5. |
| 4 | `code/assembly.py` | `assemble(problem,N,face_rule) -> A,b,grid,faces`. Interior COO triplets, then CSR. Vectorize bulk node/edge work; preserve every boundary diagonal contribution and shift known values to RHS. No dense fine-grid conversion. |
| 5 | `tests/test_structure.py` | N=3 full matrix, endpoint and central mapping, no row-wrap neighbors, nnz, symmetry, signs, positive eigenvalue on a tiny case, nonzero Dirichlet row, constant-a scaling. |
| 6 | `code/solvers.py` | Direct solver and residual-stopped GS; return solution, true residuals, success/cap flags, iterations and timers in a common result object. GS triangular split described below. |
| 7 | `code/diagnostics.py` | Named norms, truncation residual, backward error, proven spectral lower bound, estimated eigenvalues, fluxes, material-side reconstruction, bounds, rates, fit sensitivity. |
| 8 | `code/experiments.py` | E0-E4 drivers with explicit parameters and immutable prediction records; E5 optional. Require audit/check gates before emitting accepted report tables. |
| 9 | `code/plotting.py` | Read saved tables, generate surfaces, spy/error/flux/refinement/tolerance/cost figures. Label actual y coordinate for odd-N slices. |
| 10 | `run_all.py` | Extend current setup driver: execute checks, numerical experiments, metadata, CSV/JSON tables, figures and gates in that order. Add `--quick`, `--experiment E2`, `--replot`; full default must reproduce all submitted evidence. |

For GS write `A=D+L+R`, where L and R are strictly lower/upper triangular. One sweep solves `(D+L)U_new=b-R*U_old` via `scipy.sparse.linalg.spsolve_triangular`. Precompute the split outside the loop. Compute a genuine residual at the declared monitoring frequency; the last saved residual must be recomputed for the returned iterate. Return failure at iteration/time caps. No iterate-difference stopping. A zero initial residual must return success with zero sweeps and no division by zero.

The direct method delegates sparse factorization to `spsolve`; the team implements the PDE discretization and validates the delegated solve. If separating factorization and solve timings with `splu`, label the backend and ordering and cross-check against `spsolve` before grouping them as direct results.

Use a `main()` guard and package/import conventions compatible with running from the repository root. Avoid hidden working-directory assumptions. Matplotlib should use a noninteractive backend for the full reproduction command.

## Saved outputs and schemas

`data/` holds CSV/JSON numerical outputs; `figures/` holds PDF plus PNG figures. Preserve tables that support report claims in Git or in the final archive, and store the command/config/revision that regenerated them.

Minimum schemas:

- `e1_refinement.csv`: problem, N, h, dof, nnz, solver, requested_tol, flag, iterations, E2h, Einf, tau2h, tauinf, residual2rel, residualinf, backward_error, algebraic_bound2h, pairwise_order2, pairwise_orderinf, assembly_s, setup_s, solve_s.
- `e2_interface.csv`: N, face_rule, contrast, xi, E2h, Einf, crossing_flux_rel_error, side_jump_rel, left_flux_rel_error, right_flux_rel_error, minU, maxU, asymmetry, residuals, solver flag.
- `e3_tolerance.csv`: N, tol, success, sweeps, residuals, direct_difference2h, direct_differenceinf, direct_error2h, QoI, algebraic_bound2h, gate_pass.
- `e4_cost.csv`: experiment family, N, rho, solver, repeat, setup_s, solve_s, total_s, sweeps, residual_checks, requested_tol, achieved_residual, error, flux_error, nnz, sparse_bytes, factor_nnz if available, success, accuracy_gate.
- `metadata.json`: run identifier, command, git commit and dirty status, config, Python/packages, machine/thread information, timing policy, start/end times.

Suggested figures: `e0_matrix_spy`, `e1_solution_and_error`, `e1_orders`, `e2_temperature_and_flux`, `e2_error_and_bounds`, `e3_solver_floor`, `e4_cost_accuracy`, `e4_mesh_and_contrast`. Each figure has one claim, named norms, units or nondimensional labels, legend and an informative caption. Distinguish measured data from reference-slope guide lines.

## Meaningful numerical tests and falsification

These tests are appropriate because the future scientific claims depend on them; they are not yet implemented.

1. Constant-a N=3 hand matrix, boundary RHS and index map. Compare all entries, not one random action.
2. Polynomial exactness and nonzero-boundary exponential P5.
3. Forcing independently differentiated and finite-difference spot-checked; bank self-check alone is not solver verification.
4. Smooth variable-a refinement across four meshes, expected slope tested only above the algebraic floor.
5. Harmonic layered nodal/flux exactness and `[0,1]` bounds; arithmetic intentionally fails exactness while retaining symmetry and bounds.
6. Flux cancellation for shared faces and independent side-flux/exact-flux comparisons. A shared-face zero jump is not the physical-accuracy assertion.
7. Solver flag/cap behavior, initial residual zero case, true residual stopping, tighter-reference stability and error bound.
8. Falsification with missing coefficient-gradient forcing; optionally a deliberately wrong boundary sign. Keep rejected outputs distinct from accepted evidence.

For a future FEM oracle: derive the weak form and trial/test spaces, specify quadrature and element family, condense Dirichlet DoFs instead of only replacing rows, and state how the exact solution and FEM function are evaluated in a common error norm. scikit-fem availability alone is not this implementation.

## Completion sequence

1. Team reviews and owns the problem definition and E0/E1 predictions.
2. Implement stages 1-5; pass hand checks before running refinements.
3. Implement direct solver and E1/E2; create the first actual four-mesh table and interface plots.
4. Implement GS/E3 and pass algebraic-error gates; complete five-gate AI audit before using results.
5. Run E4 serially; record failed/capped cases and any absent crossover.
6. Add E5 only if E0-E4 are defensible and reproducible.
7. Write report around actual evidence, build slides from the storyboard, and have a teammate reproduce a reported number to three digits.
8. Final command regenerates all tables/figures; submit LaTeX PDF, slide PDF and runnable code zip with README, AI log and ICS.
