# MTH321 Project 2: Variable Diffusion

**Direction 3:** conservative numerical solution of `-div(a grad u) = f`.

The proposed story is: **Can a computed temperature field be accurate, physically conservative, and affordable?** Verify smooth accuracy first; expose wrong interface physics next; then compare solver cost only after algebraic error is controlled.

This repository currently contains the experiment design, theoretical predictions, implementation specification, and a verified course environment. The project experiments and final report/slides are planned work, not completed results.

- [Experiment design and theory](notes/experiment_plan.md)
- [Problem definition and conventions](notes/problem_definition.md)
- [Code implementation specification](notes/implementation_plan.md)
- [Report and presentation story](notes/report_and_slides_outline.md)
- [Course requirements mapped to evidence](notes/requirements_matrix.md)
- [Environment evidence](notes/evidence/environment_check.txt)
- [Prediction and AI audit records](notes/audit_and_predictions.md)

## Environment

Python 3.12; NumPy, SciPy, Matplotlib, ipykernel, scikit-fem, and PyAMG. The last two are required for everyone's FEM training even though this project's primary method is FDM.

Windows PowerShell, from this directory:

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.lock.txt
.\.venv\Scripts\python.exe run_all.py --setup-only
```

For an intentional dependency upgrade use `requirements.txt`, rerun the checks, and regenerate the lock with `python -m pip freeze`. The lock records the supplied environment plus the two added course packages. The `.venv` is local and excluded from Git; teammates recreate it using the lock. VS Code selects `.venv/Scripts/python.exe`.

`run_all.py` currently runs **environment checks only**, writes actual evidence, and explicitly reports that experiment generation is not implemented. Once the modules in the implementation plan are built, extend it into the brief's one-command figure/table reproduction entry point. Do not submit this setup-stage repository as a completed numerical project.

## Layout

```text
code/                  future project solver and experiment modules
  course_checks/       two unchanged supplied environment-check scripts
data/                  generated numerical tables and metadata
figures/               generated publication figures
notes/                 design, derivations, predictions, audit and evidence
report/                future LaTeX report and PDF (at least 8 pages)
slides/                future presentation source and PDF (10 minutes + Q&A)
tests/                 future matrix, boundary, flux and falsification checks
run_all.py             setup entry point; future full reproduction driver
requirements.txt       direct dependencies
requirements.lock.txt  exact installed package versions
```

Empty folders include `.gitkeep`. Course PDFs are kept outside this repository; only the two required course check scripts are copied, with provenance in `code/course_checks/README.md`. No Project 1 source or Git history is reused.

## Course actions still requiring team participation

Each member must run the checks, post their own Sanity-Check Database evidence, complete the individual challenge, and make their own Git contribution. Team members must independently audit and explain numerical code, maintain AI transparency and individual-contribution records, and reproduce reported numbers. Repository setup does not complete those personal/course submissions.
