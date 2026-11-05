# Passive-scalar channel solver

The solver uses native **parallel HDF5**, FFTW3 MPI, MPI/OpenMP, and GSL.
ESIO is no longer required. The HDF5 implementation and GSL compatibility
adapter are adapted from the neighboring `CHANNEL` project.

## Build

With GNU Fortran/C++, MPI, parallel HDF5 (including the Fortran bindings),
FFTW3 MPI, and GSL development packages installed:

```sh
make clean
make -j4
```

The default `MACHINE=LOCAL` profile uses `h5pfc` and `mpicxx`, producing
`POONGBACK_Channel_LOCAL`. Use HDF5 and FFTW built with the same MPI toolchain.
Override `H5PFC`, `FFTW_INC`, `FFTW_LIB`, `GSL_DIR`, or `GSL_LIB` for custom
installations. Historical machine profiles remain selectable with
`make MACHINE=ZHUHAI` (or another existing profile); update their dependency
paths to match the cluster. Those profiles now link `hdf5_fortran` instead
of ESIO. Always clean first: the original archive contains compiled objects,
module files, and an executable from a different installation.

## I/O and GSL compatibility

Field, statistics, spectra, and all three restart paths use
`channel_hdf5_io`. Passive-scalar datasets (`theta`, `T_mean`, `Pr`, and scalar
statistics) retain their original names. Complex fields use the historical
three-dimensional HDF5 dataspace with two-double array elements, rather than
adding a fourth scalar dimension. The reader accepts legacy ESIO layouts
0, 1, and 2 and rejects incompatible types, dimensions, and metadata.
Collective I/O supports ranks with empty selections; rank-zero statistics
use their original communicator.

`gsl_bspline_compat.h` selects the separate derivative workspace for GSL 1.x
and the merged workspace API for GSL 2.x, preserving the numerical calls.
See the [GSL 2.0 release notes](https://lists.gnu.org/archive/html/info-gnu/2015-10/msg00014.html)
for the API change.

## Validation

Run on a Linux MPI installation with Python, NumPy, and h5py:

```sh
python3 tests/run_build_tests.py
python3 tests/io/run_io_tests.py
```

The build test performs clean debug and optimized parallel builds and checks
B-spline partition of unity, derivative sums, and integration on uniform
and stretched grids. The I/O test checks exact field values including theta,
mean scalar profiles, complex datatypes, spectra, append behavior, different
writer/reader decompositions, empty ranks, and malformed input rejection.
Legacy fixtures were generated with real ESIO and copied from `CHANNEL/tests/io/fixtures`.
These tests validate the dependency migration; they do not establish a full
DNS trajectory regression for this passive-scalar solver.
