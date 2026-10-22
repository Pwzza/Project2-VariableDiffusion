# Requirements, sources and planned evidence

Local sources reviewed on 2026-10-22; source pages below refer to PDF page numbering. Paths are relative to the parent of this repository unless otherwise stated.

| Source | Requirement/inspiration | Implementation or planned evidence |
|---|---|---|
| `Materials/01_Project_Brief/Project_2_Brief.pdf`, pp.1-4 | Conservative FDM or justified FEM; explicit boundary treatment; four meshes and named norms; solver-error separation; substantial depth; reproducibility | Definition sheet; E0-E4; implementation plan; future full `run_all.py` |
| Same, p.2 | >=8-page LaTeX report PDF; code zip/README; slide PDF, 10 min + 3 min Q&A | Report/storyboard outline; report/ and slides/ reserved |
| Same, pp.2,4 | At least two independent checks; predictions before solving; AI audit and transparency; ICS | Analytic smooth/layered references, hand derivations, direct-vs-GS; audit record and future team verification |
| `Materials/01_Project_Brief/problem_pack.pdf`, Direction 3, pp.5-7 | Prescribed smooth MMS, independently confirmed forcing, four-mesh L2/Linf orders and residuals | E1, manual f derivation, independent finite-difference check |
| Same, pp.5-7 | Shared positive faces, correct resistance, Dirichlet elimination, asymmetry/SPD evidence | Definition sheet energy proof; E0 checks, eigenvalue corroboration |
| Same, pp.5-7 | Layered arithmetic/harmonic comparison, odd interval grids, temperature and flux errors/jump, maximum-principle check | E2; independent material-side flux and exact-flux comparisons; bounds |
| Same, p.7 | Conditioning or iterations under mesh refinement and contrast at comparable residual tolerance | E4 with verified solution/flux accuracy as an additional fairness gate |
| Same, Section 4, p.7 | One substantial matched comparison with stated cost measure | Choose direct vs GS (item 2), in addition to core interface comparison; consider other items explicitly |
| Same, Section 5, pp.8-10 | Optional original extensions | Defer L-shape/adaptivity/Neumann/nonlinear scope; E5 is a related modest extension |
| `Materials/04_Tutorials/environment_setup_bvp.pdf`, pp.1-2 | All six packages including skfem/pyamg; supplied complete setup check; bank self-check; evidence | Existing .venv completed, frozen dependencies, actual `notes/evidence/environment_check.txt` |
| `tutorial5_slides.pdf`, pp.7-13 | Individual P5 reference challenge; hand matrix; Git bring-up; code/notes/report layout; personal course records | E0 planned; repo layout created; personal submissions remain team actions |
| `fdm_worksheet.pdf`, pp.1-3 | Signed stencil, exact N=3 matrix, boundary shifts, sparse assembly, separate residual/truncation/error | Definition sheet plus E0/E1/E3 |
| `tutorial6_slides.pdf`, `scipy_sparse_cheatsheet.pdf` | Sparse storage, convergence fits, GS/CG/AMG and scaling | Implementation plan and E4; CG/AMG optional additions |
| `mms_worksheet.pdf`, pp.1-2; `tutorial7_slides.pdf` | Finest fitted mesh algebraic gate <=1e-3 of discretization error; tighter reference; at least 3 asymptotic levels; falsification | E1/E3; use brief's stronger direct-error requirement even though worksheet allows just reporting direct residual |
| `precomputation_questions_worksheet.pdf`, pp.1-2 | Predict -> run -> reconcile, quantify signs/symmetries/limits and error sources | Immutable pre-run predictions plus append-only outcomes in audit document |
| `verification_clinic.pdf`, `session_map.pdf` | PDE, representative row, quantitative result, tolerance and unresolved issue; readiness gates | Definition sheet/E0; future four-mesh result; audit status must be honest |
| `ai_solver_audit_worksheet.pdf`, pp.1-2 | Five gates before numerical solver output is used | `audit_and_predictions.md`; do not claim team review has happened |
| `tutorial8_slides.pdf`, `feedback_to_evidence_worksheet.pdf` | Convert feedback into quantified changes; limits and cost | Feedback ledger and E4; record table/figure/number that resolves each comment |
| `tutorial9_slides.pdf` | Reproduce a reported number to three digits, reconcile discrepancies, final package | Reproduction ledger; future teammate clean-environment run |
| `fem_worksheet.pdf`, `nonrectangular_domains_guide.pdf` | Weak-form/mesh/space/quadrature/boundary clarity if using FEM or irregular domains | FEM alternative scoped in plan; not silently claimed as implemented |
| `../Project 1 - Completed/NumericalPDE/report/slide.pdf`, slides 2-4,6-29 | Story-driven staging: sanity, exact reference, structural stress, reference reliability, cost | E0-E4 narrative and 13-slide storyboard; no SDE conclusions transferred to elliptic PDEs |

Missing in the provided tree: the brief's `06_Templates` and the shared `git_quickref` guide. The existing Project 1 repository and Tutorial 5 establish the familiar code/notes/report layout; data/figures/slides/tests support this project's reproduction plan. `.gitkeep` preserves empty folders. No missing document is treated as read.

## Required scope versus optional scope

Core: conservative Direction 3 scheme, smooth refinement, interface mean comparison, boundary bounds, symmetry/SPD, mesh/contrast solver evidence. Selected depth: direct vs GS at matched accuracy with stated cost. Optional: CG/AMG, off-center interface E5, nonzero-boundary variable-a MMS, or a fully matched 2D FEM method comparison. No need to complete every depth/open direction.

## Course milestone actions

Week 5: problem definition Parts 1 and 2 for FDM; hand matrix; per-person sanity database and individual challenge; team topic/weekly entry and forum record; every member makes a Git contribution.

Week 6: working code, first MMS refinement table, chosen depth comparison.

Week 7: clinic readiness sheet with four-mesh table, named norm, tolerances, provisional slope and remaining question.

Week 8: defensible comparison and limitations; if volunteering for first presentation wave, deliver work-in-progress slides and draft.

Week 9/common deadline: second wave if selected; final report with AI log/ICS, slide PDF and code archive. The brief states Monday of Week 10 at 23:59. Determine the calendar date from the course calendar; do not infer it from today's date.
