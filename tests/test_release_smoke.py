from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def _clean_runtime_caches():
    """Remove caches created by the test runner itself before repository audit."""
    shutil.rmtree(ROOT / '.pytest_cache', ignore_errors=True)
    for p in ROOT.rglob('__pycache__'):
        shutil.rmtree(p, ignore_errors=True)
    for p in ROOT.rglob('*.pyc'):
        try:
            p.unlink()
        except FileNotFoundError:
            pass


def test_repo_doctor():
    _clean_runtime_caches()
    subprocess.run(
        [sys.executable, str(ROOT / 'scripts/repo_doctor.py'), '--repo', str(ROOT), '--strict'],
        check=True,
    )


def test_checksums():
    subprocess.run(
        [sys.executable, str(ROOT / 'scripts/verify_release.py'), '--repo', str(ROOT)],
        check=True,
    )
