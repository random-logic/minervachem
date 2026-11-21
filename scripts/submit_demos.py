"""Run the shared submitter and report Slurm rejection details to the user."""
from pathlib import Path
import runpy
import subprocess


if __name__ == '__main__':
    runner = Path(__file__).resolve().parents[2] / 'buhito/src/buhito/demo_runner.py'
    try:
        runpy.run_path(str(runner), run_name='__main__')
    except subprocess.CalledProcessError as error:
        detail = (error.stderr or error.stdout or '').strip()
        raise SystemExit(
            f'Slurm submission failed (exit {error.returncode}):\n'
            f'{detail or "sbatch returned no diagnostic"}'
        ) from None
