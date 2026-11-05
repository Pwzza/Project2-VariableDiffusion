from pathlib import Path
import argparse
import contextlib
import importlib
import io
import platform
import sys
import os
import subprocess

ROOT = Path(__file__).resolve().parent


def bootstrap_environment():
    
    if sys.prefix != sys.base_prefix:
        return

    venv_dir = ROOT / ".venv"
    if os.name == "nt":  # Windows
        python_exe = venv_dir / "Scripts" / "python.exe"
        pip_exe = venv_dir / "Scripts" / "pip.exe"
    else:                # macOS / Linux
        python_exe = venv_dir / "bin" / "python"
        pip_exe = venv_dir / "bin" / "pip"

    if not venv_dir.exists():
        print("Setting up virtual environment...")
        import venv
        venv.create(venv_dir, with_pip=True)

    req_file = ROOT / "requirements.lock.txt"
    if not req_file.exists():
        req_file = ROOT / "requirements.txt"

    if req_file.exists():
        print(f"Installing dependencies from {req_file.name}...")
        subprocess.check_call([str(pip_exe), "install", "--upgrade", "pip", "-q"])
        subprocess.check_call([str(pip_exe), "install", "-r", str(req_file)])
    else:
        print("Warning: No requirements file found.")

    print("Environment ready. Restarting script inside the virtual environment...\n")
    os.execv(str(python_exe), [str(python_exe)] + sys.argv)

bootstrap_environment()
sys.path.insert(0, str(ROOT / 'code'))


def run_script(path):
    path = Path(path)
    old_path = sys.path.copy()
    sys.path.insert(0, str(path.parent))
    try:
        with open(path, 'r', encoding='utf-8') as f:
            code = compile(f.read(), str(path), 'exec')
        namespace = {'__name__': '__main__', '__file__': str(path)}
        exec(code, namespace)
    finally:
        sys.path[:] = old_path


def setup_checks():
    log = io.StringIO()
    with contextlib.redirect_stdout(log):
        print('MTH321 Project 2 environment evidence')
        print('Python:', sys.version.replace('\n', ' '))
        print('Platform:', platform.platform())
        for name in ('numpy', 'scipy', 'matplotlib', 'ipykernel', 'skfem', 'pyamg'):
            module = importlib.import_module(name)
            print(f'{name}: {module.__version__}')
        for filename in ('poisson_fd_setup_check.py', 'test_problem_bank.py'):
            print(f'\nCourse check: {filename}')
            run_script(ROOT / 'code' / 'course_checks' / filename)
        print('\nBank result: 0 failures (supplied self-check returned without assertion failure).')
        print('Setup checks: PASS. This does not verify the future variable-diffusion solver.')
    output = log.getvalue()
    path = ROOT / 'notes' / 'evidence' / 'environment_check.txt'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(output, encoding='utf-8')
    print(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--setup-only', action='store_true', help='run the supplied environment checks')
    parser.add_argument('--quick', action='store_true', help='four smaller meshes for development; not final evidence')
    parser.add_argument('--plots-only', action='store_true', help='regenerate figures from existing data')
    parser.add_argument('--cost-only', action='store_true', help='benchmark saved cases and redraw figures')
    args = parser.parse_args()
    
    if args.setup_only:
        setup_checks()
    else:
        from diffusion.experiments import run
        from diffusion.plots import generate
        if not args.plots_only and not args.cost_only:
            run(ROOT, quick=args.quick)
        if not args.plots_only:
            from diffusion.presentation_studies import benchmark
            benchmark(ROOT)
        generate(ROOT)
        print('Done: numerical data in data/ and presentation figures in figures/.')