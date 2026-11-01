# MTH321 Project 2: Variable Diffusion

**Direction 3:** conservative numerical solution of `-div(a grad u) = f`.

The proposed story is: **Can a computed temperature field be accurate, physically conservative, and affordable?** Verify smooth accuracy first; expose wrong interface physics next; then compare solver cost only after algebraic error is controlled.

The numerical code now implements **both conservative FDM and two-dimensional P1 FEM**, four-mesh smooth/interface experiments, diagnostics, and reproducible presentation figures. The final report and Beamer slides are not assembled yet.

**Start with [the current code guide](notes/code_guide.md).** The October planning documents linked below preserve the original design; their standalone sanity/error/cost experiment numbering and FDM-only emphasis have been superseded by the approved two-model, two-method presentation plan.

- [Experiment design and theory](notes/experiment_plan.md)
- [Problem definition and conventions](notes/problem_definition.md)
- [Code implementation specification](notes/implementation_plan.md)
- [Report and presentation story](notes/report_and_slides_outline.md)
- [Course requirements mapped to evidence](notes/requirements_matrix.md)
- [Environment evidence](notes/evidence/environment_check.txt)
- [Prediction and AI audit records](notes/audit_and_predictions.md)

## Environment

Python 3.12; NumPy, SciPy, Matplotlib, ipykernel, scikit-fem, and PyAMG. scikit-fem implements the project's 2D P1 FEM pathway. PyAMG is available for optional future work; the current iteration comparison uses unpreconditioned CG.

Windows PowerShell, from this directory:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.lock.txt
.\.venv\Scripts\python.exe run_all.py --setup-only
```

For an intentional dependency upgrade use `requirements.txt`, rerun the checks, and regenerate the lock with `python -m pip freeze`. The lock records the supplied environment plus the two added course packages. The `.venv` is local and excluded from Git; teammates recreate it using the lock. VS Code selects `.venv/Scripts/python.exe`.

Run from this directory (no activation needed):

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -v
.\.venv\Scripts\python.exe run_all.py
```

The default command runs the full four-mesh experiments, checks acceptance criteria, and writes JSON/CSV/NPZ data plus PDF/PNG figures. `--quick` uses smaller four-mesh studies for development. It overwrites the same output names: regenerate with the default command before using final figures. `--plots-only` redraws figures from saved data. `--setup-only` retains the original course environment checks.

## Layout

```text
code/diffusion/        models, sparse solvers, diagnostics, experiments, plots
  course_checks/       two unchanged supplied environment-check scripts
data/                  generated tables, fields, profiles, checks and metadata
figures/               generated publication figures
notes/                 design, derivations, predictions, audit and evidence
report/                future LaTeX report and PDF (at least 8 pages)
slides/                future English Beamer source/PDF (target 8 min 30 sec)
tests/                 independent hand-matrix, boundary, flux and convergence checks
run_all.py             full experiment and figure reproduction driver
requirements.txt       direct dependencies
requirements.lock.txt  exact installed package versions
```

Empty folders include `.gitkeep`. Course PDFs are kept outside this repository; only the two required course check scripts are copied, with provenance in `code/course_checks/README.md`. No Project 1 source or Git history is reused.

## Course actions still requiring team participation

Each member must run the checks, post their own Sanity-Check Database evidence, complete the individual challenge, and make their own Git contribution. Team members must independently audit and explain numerical code, maintain AI transparency and individual-contribution records, and reproduce reported numbers. Repository setup does not complete those personal/course submissions.
