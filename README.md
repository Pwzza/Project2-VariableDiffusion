# MTH321 Project 2: Variable Diffusion - Team 11

**Conservative numerical solution of `-div(a grad u) = f`.**

The following project attempts to answer: **Can a computed temperature field be accurate, physically conservative, and affordable?** 

The numerical code implements **both conservative FDM and two-dimensional P1 FEM**, four-mesh smooth/interface experiments, diagnostics, cost studies, and reproducible presentation figures. The expanded Beamer presentation is available in `slides/main.tex` and `slides/main.pdf`; the final report remains pending.

Our team automated eveything, so the script will activate all the environments and install the necessary packages. 

### Linux / MacOS
```bash
cd Project2-VariableDiffusion && python3 run_all.py
```

### Windows
```powershell
cd Project2-VariableDiffusion && python run_all.py
```
