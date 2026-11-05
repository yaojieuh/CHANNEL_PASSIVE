"""Clean builds and GSL numerical checks, independent of archived object files."""
from pathlib import Path
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def run(args, cwd):
    result = subprocess.run(args, cwd=cwd, text=True, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, timeout=180)
    if result.returncode:
        raise RuntimeError(result.stdout)
    return result.stdout


def main():
    with tempfile.TemporaryDirectory(prefix='passive-build-') as temp:
        folder = Path(temp)
        for path in ROOT.iterdir():
            if path.suffix in ('.f90', '.cpp', '.h') or path.name == 'makefile':
                shutil.copy2(path, folder/path.name)
        for mode, flags in [('debug', '-O0 -g -fopenmp -fallow-argument-mismatch -fcheck=all'),
                            ('optimized', '-O3 -fopenmp -fallow-argument-mismatch')]:
            run(['make', 'clean'], folder)
            run(['make', '-j4', f'FFLAGS={flags}'], folder)
            print(f'PASS: clean parallel {mode} solver build', flush=True)
        run(['g++', '-I'+str(ROOT), str(ROOT/'labs_fill_full.cpp'),
             str(ROOT/'tests/test_bspline.cpp'), '-lgsl', '-lgslcblas',
             '-o', 'test_bspline'], folder)
        print(run([str(folder/'test_bspline')], folder), end='')


if __name__ == '__main__':
    main()
