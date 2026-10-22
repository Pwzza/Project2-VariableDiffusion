"""Setup-stage entry point. No project experiment results exist yet."""
from pathlib import Path
import argparse
import contextlib
import importlib
import io
import platform
import runpy
import sys

ROOT = Path(__file__).resolve().parent


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
    args = parser.parse_args()
    setup_checks()
    if not args.setup_only:
        print('PLANNED: E0-E5 figure/table generation is not implemented. See notes/implementation_plan.md.')
