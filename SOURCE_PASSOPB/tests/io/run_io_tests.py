"""Verify values, schema, empty selections, legacy layouts, and rejection."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]


def run(command, cwd, expected=None):
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True, timeout=90)
    output = result.stdout+result.stderr
    if expected:
        assert result.returncode != 0 and expected in output, output
    else:
        assert result.returncode == 0, output


def main():
    launcher = ['mpiexec', '--allow-run-as-root', '--oversubscribe']
    cases = 0
    with tempfile.TemporaryDirectory(prefix='channel-io-') as temp:
        folder = Path(temp)
        exe = folder/'io_roundtrip'
        run(['h5pfc', '-O0', '-g', '-fcheck=all', str(ROOT/'module_hdf5_io.f90'),
             str(ROOT/'tests/io/io_roundtrip.f90'), '-o', str(exe)], folder)
        for writers in (1, 2, 4):
            path = folder/f'native-{writers}.h5'
            run([*launcher, '-n', str(writers), str(exe), 'write', str(path)], folder)
            cases += 1
            with h5py.File(path) as h:
                for field, offset in [('u', 1), ('v', 2), ('w', 3), ('g', 4), ('phi', 5), ('theta', 6)]:
                    assert h[field].shape == (4, 5, 3)
                    assert h[field].dtype.subdtype == (np.dtype('float64'), (2,))
                    z, x, y, c = np.indices((4, 5, 3, 2))
                    want = (c+1)+10*(y+1)+100*(x+1)+1000*(z+1)+10000*offset
                    np.testing.assert_array_equal(h[field][()], want)
                assert h['appended'][0] == 43
                np.testing.assert_array_equal(h['Euu_kx'][()],
                    (np.arange(3)+1)[None, :] + 100*(np.arange(5)+1)[:, None])
            for readers in (1, 2, 4, 6):
                run([*launcher, '-n', str(readers), str(exe), 'read', str(path)], folder)
                cases += 1
            print(f'PASS: {writers} writers, 1/2/4/6 readers; exact values, array dtype and empty selections', flush=True)
        for fixture in sorted((ROOT/'tests/io/fixtures').glob('*.h5')):
            for readers in (1, 2, 4, 6):
                run([*launcher, '-n', str(readers), str(exe), 'read', str(fixture), 'legacy'], folder)
                cases += 1
            print(f'PASS: real legacy fixture {fixture.name}, 1/2/4/6 readers', flush=True)
        # Reject ambiguous layouts/types instead of treating arbitrary bytes as complex data.
        for case in ('rank', 'type', 'shape', 'metadata'):
            path = folder/f'bad-{case}.h5'
            shutil.copyfile(folder/'native-1.h5', path)
            with h5py.File(path, 'r+') as h:
                data = h['u'][()]
                if case == 'metadata':
                    h['u'].attrs['esio_field_metadata'] = np.zeros(7, dtype=np.int32)
                else:
                    del h['u']
                    if case == 'rank':
                        h['u'] = data
                    elif case == 'type':
                        h['u'] = data[..., 0]
                    else:
                        h.create_dataset('u', shape=(4, 6, 3), dtype=np.dtype((np.float64, (2,))))
            expected = {'rank': 'unexpected dataset rank', 'type': 'expected array-valued field',
                        'shape': 'restart dimensions do not match input', 'metadata': 'malformed legacy metadata length'}[case]
            run([*launcher, '-n', '4', str(exe), 'read', str(path)], folder, expected)
            cases += 1
        print(f'PASS: {cases} MPI I/O executions, including four expected rejections', flush=True)


if __name__ == '__main__':
    main()
