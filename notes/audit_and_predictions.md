# Prediction, audit and reproduction records

## Current status

This setup-stage repository includes AI-assisted planning and an environment runner. The two course numerical checks are copied unchanged with attribution. No project variable-diffusion solver or claimed experimental result has yet been generated. The design's mathematical predictions are reasoned derivations for team review; they do not replace independent ownership or validation.

## Predictions written before project runs

| ID | Prediction | Basis | Output/reconciliation |
|---|---|---|---|
| P0 | Constant-a N=3 matrix equals the stated 4x4 matrix; polynomial reproduced up to round-off | Central differences on a quadratic; boundary elimination | Pending project implementation |
| P1 | Smooth E2h/Einf and scaled truncation residual approach h^2 | Taylor consistency plus stability for smooth uniformly elliptic square problem | Pending |
| P2 | Harmonic layered U is exact at nodes and q=-20/11, up to solver error | Exact subedge resistance and piecewise linear profile | Pending |
| P3 | Arithmetic layered closure can remain SPD/bounded yet be inaccurate | Positive shared weights do not encode the right resistance | Pending |
| P4 | Missing coefficient-gradient forcing produces a persistent error against intended u* | It changes the differential equation | Pending falsification |
| P5 | Tightening solver tolerance eventually leaves PDE error unchanged | Algebraic/discretization decomposition | Pending |
| P6 | Refinement generally increases stationary solver work; contrast can worsen conditioning | Spectral comparison and stationary iteration analysis | Pending; no monotonic contrast or crossover guarantee |
| P7 | Exact distance-weighted resistance reproduces an off-center piecewise-linear benchmark | Integrating du=C dx/a on the crossing edge | Optional, pending |

Preserve this table's initial predictions. Append measured outputs and reconciliation decisions; do not retrospectively rewrite a prediction to match data.

## Five-gate project solver audit

| Gate | Required evidence | Status |
|---|---|---|
| 1 Transcription | PDE sign, a/f/g/interface parameters and stencil checked term by term | Planned; team review pending |
| 2 Structure/indexing | Hand rows, N=3 matrix, boundary/interface rows, ij index map, sparse structure | Pending solver implementation |
| 3 Linearisation | Nonlinear Jacobian test | Not applicable to this linear PDE; stopping/termination still audited in E3 |
| 4 Independent reproduction | Exact references, predicted refinement order, independent solver or hand limit; provenance of each result | Pending experiments |
| 5 Interrogation | Attack hypotheses and tests below, followed by team sign-off | Hypotheses recorded; tests pending |

Three ways plausible convergence can conceal a defect:

1. Converging to the wrong operator because f omits grad(a) dot grad(u): independently derive f and run E1's falsification.
2. Symmetric bounded arithmetic-interface solution with incorrect heat flow: use exact-flux and independent side-reconstruction checks in E2.
3. Tiny normalized residual hiding forward error from conditioning or large boundary b: use E3 residual-to-error bound and tighter independent solve, especially on the finest mesh/large contrast.

Other threats: transposed plot axes, wrong all-boundary g, shared-face zero-jump tautology, residual-monitoring cost bias, a solve silently hitting its iteration cap, fitting a line through round-off data. Each has a specific check in the implementation/experiment plan.

## AI Transparency Log

| Date | AI contribution | Human change/review | Validation | Status |
|---|---|---|---|---|
| 2026-10-22 | Read supplied materials; drafted experiment design, theory, problem definition, implementation and presentation outline; completed environment folders/dependencies and setup runner | To be completed by team; do not invent reviewer names or contributions | Supplied setup scripts executed; see environment evidence. Scientific predictions and project solver still require independent review/experiments | Setup only |

## Feedback-to-evidence ledger

| Observation | Claim affected | Action | Numerical evidence / location | Remaining risk |
|---|---|---|---|---|
| Pending clinic/team feedback | | | | |

## Independent reproduction ledger

| Report number | Generating command/config/revision | Reproduced output | Match to 3 digits? | Reviewer and date |
|---|---|---|---|---|
| No report numbers yet | | | | |

Do not prefill team sign-off, individual contribution statements, course database submissions or reproduction successes. Each member must run, review, explain and record their actual work.
