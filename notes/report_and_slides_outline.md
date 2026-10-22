# Report and presentation design

This is an outline for future final deliverables. It follows Project 1's motivation -> prediction/theory -> implementation -> evidence rhythm. Predicted behavior must remain labelled as a prediction until actual runs and audit support it. Use the previous deck's staged story structure; do not repeat its 32-slide length for a 10-minute talk.

## Report: planned 10-12 main-text pages plus appendices

The brief requires **at least 8 pages in LaTeX, submitted as PDF**. The referenced `Materials/06_Templates/latex_report_template_bvp.tex` is absent from the supplied materials; use the specified section order and add the original course template if later obtained.

| Section | Approx. pages | Purpose and quantitative evidence |
|---|---:|---|
| Background | 1 | Heat transport through layered composites; why temperature, flux, and computational reliability all matter. State the question and scope. |
| Model | 1 | PDE, domain, bounds on a, sign convention, perfect-interface conditions, all-boundary Dirichlet data; exact smooth and layered references. |
| Discretisation | 2 | Flux derivation, resistance average, COO/CSR assembly, index map, Dirichlet elimination, tiny matrix, symmetry/SPD/M-matrix proof. |
| Verification | 2 | E0, independent forcing check, named E1 refinement tables, residual/truncation/error separation, E3 algebraic gate and falsification. |
| Numerical experiments | 3-4 | E2 temperature/flux comparison and bounds; E4 matched-accuracy cost, mesh/contrast conditioning and iterations; E5 if completed. Report failed expectations and limits. |
| Conclusions | 0.5-1 | Only claims backed by evidence; practical method/solver choice under tested conditions, unresolved issues and possible extensions. |

Appendices: full derivations, full tables and parameters, five-gate code audit, pre-computation reconciliations, AI Transparency Log, Individual Contribution Statements, and reproduction ledger. Page allocations are a writing guide, not a reason to pad the report.

Every main result must name the generating command/config, data table, reference, error norm, tolerance and verification gate. Figures should interpret numbers rather than substitute for them. Clearly distinguish mandatory core, selected depth and optional extensions.

## Presentation: 13 main slides, 10 minutes + 3 minutes Q&A

Each title below is the proposed **takeaway**, followed by the visual/evidence and transition. Replace conditional/predicted wording with measured claims only after results exist. Time budget totals 600 seconds.

| # | Time | Takeaway title | Content and spoken connection |
|---|---:|---|---|
| 1 | 25 s | A temperature field is useful only if its heat flow can be trusted | Title, team, composite sketch, practical question. |
| 2 | 40 s | Reliability requires accuracy, correct interface transport, and controlled cost | E0-E4 story map, matching Project 1's staged experimental structure. Explain why the checks build on each other. |
| 3 | 45 s | Divergence form expresses conservation across a material interface | PDE, q=-a grad u, model assumptions and interface conditions; schematic of continuous u and kinked slope. |
| 4 | 50 s | One shared face flux gives symmetry; resistance determines its physical value | Five-point row, boundary elimination, tiny matrix inset; derive harmonic mean from two subedge resistances. |
| 5 | 35 s | Exact simple cases catch indexing and boundary errors cheaply | E0 hand matrix, nonzero-g check and polynomial exactness; distinguish setup success from solver verification. |
| 6 | 50 s | Smooth manufactured data predict second-order temperature accuracy | E1 f derivation with gradient-coefficient term highlighted; L2h/Linf definitions and h^2 expectation. |
| 7 | 50 s | The refinement table tests the prediction, not the appearance of the plot | Actual E1 table/log-log fit, pairwise rates, residual and selected fit levels. State any deviations. |
| 8 | 50 s | A material jump is the stress test that smooth convergence cannot replace | E2 exact temperature/flux derivation; predicted aH=20/11 versus aA=11/2, matched interface setup. |
| 9 | 60 s | A bounded symmetric solution can still have the wrong heat flux | Actual E2 temperature and flux slice, numerical errors and bounds. Shared-face cancellation is not proof of correct flux magnitude. |
| 10 | 50 s | We must finish the linear solve before measuring PDE accuracy | E3 error decomposition and tolerance-floor plot, quantified algebraic-error gate. Connect to the error-source discipline of Project 1. |
| 11 | 55 s | Fair solver comparisons use matched accuracy and include setup cost | E4 timing protocol, error-vs-cost and iteration-vs-mesh/contrast plots; display failed/capped cases. |
| 12 | 50 s | The evidence determines the practical solver choice on the tested range | Measured tradeoff, contrast conditioning, crossover only if observed. Optional E5 inset if it adds a distinct finding. |
| 13 | 40 s | Reliable transport needs verified discretisation, interface physics, and solver accuracy | Three supported findings, scope limits, reproducibility command and link/QR to repo. End with conclusions that answer slide 1. |

Design: use a restrained blue/white mathematical style similar to the prior deck, large labelled plots, one main message per slide, short equations, and a small motivation/prediction/evidence marker if helpful. Do not cram four unrelated plots into one panel. Reserve derivations and full data for backup.

## Backup slides for Q&A

- Full SPD/M-matrix derivation and numerical asymmetry/eigenvalue corroboration.
- Interior/boundary/interface rows and ij vector indexing.
- Exact resistance and general contrast/interface formulas.
- Why odd-interval meshes do not halve h exactly and have no y=1/2 grid row.
- Why smooth forcing changes sign and the harmonic maximum-principle bound cannot be applied there.
- Algebraic residual-to-error bound and direct reference reliability.
- Why arithmetic full-Dirichlet solutions need not stay one-dimensional.
- Complete solver stopping/cap policy, conditioning bounds and timing scope.
- AI audit, independent reproduction, contribution record and limitations.

## Report versus slides

The report is organized by mathematical and evidential completeness. The slides follow the audience's questions: why heat flow -> how conservation is represented -> why exact tests -> why an interface stress test -> why solver error -> what accuracy costs. This preserves the useful structure of the previous presentation without treating a list of assigned tasks as a scientific motivation.
