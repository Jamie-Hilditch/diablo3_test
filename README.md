# Diablo 3

This version is based on the original version of J. R. Taylor (cf. DIABLO2.0 ReadMe).
It includes updated FFTW and pHDF5 libraries. It also requires the _toml-f_ library. The appropriate paths need to be set in the _Makefile_ (_for/Makefile0_) and 
added to the `LD_LIBRARY_PATH`.

## Test case
The test case is a 10 minute simulation on a fairly small domain. The initial conditions are pretty much nonsense. 

The number of grid points in x,y,z are set in _grid_def.all_. The number of grid points in y must match _grid.h5_.

The number of processes used in the y and z directions are passed as command line arguments to _for/setup_run_. The number of processes in y must divide NY - 1 = 193 - 1 = 192 = 6 x 32 and the number of processes in z must divide NZ = 32.

In _input.toml_ the VERBOSITY is set to its highest setting. This includes printing out a message every time `mpi_alltoall` is called in `fft_xz_to_fourier`. 

## Running the test case:
1. Copy _test_case_ directory to $SCRATCH
2. Modify the path in _diablo.slurm_ that points to _for/setup_run_
3. Submit slurm job

## Important files
- _for/Makefile0_ is the makefile. We usually recompile the code every time we run it as the grid sizes and processor distribution need to be known at compile time.
- _for/setup_run_ creates the grid files included in the fortran code and then compiles the code. This also sets a number of flags that set options in the makefile. The `--toml` flag is necessary whereas the `--debug` flag is optional. We don't use the `--shared` flag. This was an attempt to improve the performance of the code using mpi shared memory directives but it turns out to be slower. Turning this option on would actually bypass the `mpi_alltoall` call that keeps failing but I tried this and this code also fails. 
- _for/domain.f90_ this is where `mpi_init`is run and the mpi_subcarts are created.
- _for/fft.f90_ this file isolates all the calls to fftw. Contains the subroutine `fft_xz_to_fourier` in which `mpi_alltoall` keeps failing. Note that `fft_xz_to_physical` also has a call to `mpi_alltoall` which doesn't seem to fail.
