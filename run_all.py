"""One-command numerical experiments and presentation figure reproduction."""
from pathlib import Path
import argparse
import contextlib
import importlib
import io
import platform
import runpy
import sys

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'code'))


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
            runpy.run_path(str(ROOT / 'code' / 'course_checks' / filename), run_name='__main__')
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
