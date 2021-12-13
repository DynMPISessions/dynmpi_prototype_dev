# Dynamic Resource Management using MPI Sessions on p4est #

This project introduces new functionality to the p4est library for parallel AMR
([www.p4est.org](https://www.p4est.org)) by C. Burstedde et al. to support adapting 
the resources used during runtime.

The changes to `p4est` are managed in their own repository, which can be found 
in the `p4est` submodule.

For dynamic resource management, the `libmpidynres` simulation layer for 
MPI Sessions written by Jan Fecht is used. For this, it was extended by ad-hoc functions 
for requesting user-defined resource changes. The code can be found in the 
`libmpidynres` submodule.

The main purpose of this repo is to provide example applications for demonstrating the 
use of the added functionality. The applications are located in the `applications/` 
direcory and are divided into three categories, the unit tests, examples for demonstrating 
the functionality, and benchmarks for evaluating the performance.

The example simulations calculate the shallow-water equations. For this, they use the
SWE solvers from the Chair of Scientific Computing at TUM, which can be found in the
`applications/swe_solvers` submodule.

Animations from different simulation runs can be found 
on [YouTube]( https://www.youtube.com/playlist?list=PLmF0b5uzp4lzxNKeycyHIOEhc7T-RWT2O)


## Compiling and Installing ##

### p4est ###

The `p4est` library can be built using the `Makefile` in `p4est/`.

The full build instructions can be found in the [p4est howto](https://p4est.github.io/p4est-howto.pdf)
or the `p4est/README.md`.
Normally, you might want to compile it by calling 

`./configure --enable-mpi && make && make install`

from the `p4est/` directory. When compiling for the first time, it might be necessary to call
`./bootstrap` in `p4est/` first.

In order to make the library accessible to the applications, the `LD_LIBRARY_PATH` has to include
the `p4est/local/lib/` directory.

### libmpidynres ###

For compiling the MPI Sessions simulation layer, simply call `make &&make install` from the
`libmpidynres/` directory. The full build options can be found in the `libmpidynres/README.md`.

If required, the build path can be altered by editing the makefile.

### Example Applications ###

The example applications use the `SCons` build system. Therefore, `SCons` has to be installed first.

After that, the code can be compiled by calling `scons [OPTIONS]` from the `applications/` directory.
The available options are:

 * `-h`           ~ Print information about the compile options
 * `buildDir`     ~ Optional parameter for the build directory, default is `applications/build/`
 * `example`      ~ Which application to build. Default value is `mpidynres2d`, possible values are:
   * `incProcs`               ~ The example simulation using the naive approach for Increasing
                                #Processes on a two dimensional adaptive grid
   * `mpidynres2d`            ~ The two dimensional example simulation using libmpidynres
                                for adaptive resource management on an adaptive grid
   * `mpidynres3d`            ~ Same as above, but in three dimensions
   * `benchConstProcsFixed`   ~ The benchmark using a Constant #Processes on a Fixed Grid
   * `benchConstProcsSplit`   ~ The benchmark using a Constant #Processes on a Split by Four Grid
   * `benchConstProcsAdapt`   ~ The benchmark using a Constant #Processes on an Adaptive Grid
   * `enchIncProcsFixed`      ~ The benchmark using an Increasing #Processes on a Fixed Grid
   * `benchIncProcsSplit`     ~ The benchmark using an Increasing #Processes on a Split by Four Grid
   * `benchMpidynresFixed`    ~ The benchmark using libmpidynres on a Fixed Grid
   * `benchMpidynresSplit`    ~ The benchmark using libmpidynres on a Split by Four Grid
   * `benchDecMpidynresFixed` ~ The benchmark using Decreasing libmpidynres on a Fixed Grid
   * `benchAdapt`             ~ The benchmark using an Adaptive #Processes on an Adaptive Grid
   * `tests`                  ~ The test suite
 * `compileMode`  ~ How to compile the examples, can be `debug` or `release`, default value is `release`
 * `netCDF`       ~ Whether to use the `netCDF` library which is required for reading scenarios
                    from files. Default value is `false`.
 * `netCDFDir`    ~ Location of the `netCDF` library, if it should be used and is not installed
                    at the default path

For example, the three-dimensional example application could be compiled in debug mode using

`scons example=mpidynres3d compileMode=debug`

Unless specified otherwise, the applicaitons will be built to the `applications/build/` directory.

The names of the executables are constructed as `SWE_p4est_[EXAMPLE_NAME]_[COMPILE_MODE]`, where
the `EXAMPLE_NAME` is the value passed to the `scons example` parameter and the `COMPILE_MODE` is
the value specified with `scons compileMode`.


## Running the Code ##


### Test Suite ###

The test suite can be executed by calling

`applications/build/SWE_p4est_tests_release`

If compiled in debug mode, the name of the executable is `SWE_pd4est_tests_debug`. 

This program does not take any options, and does not require starting it with `mpirun`.
All necessary processes are spawned using `MPI_Comm_spawn` by the appliation.


### Benchmarks ###

The benchmarks accept the following command line options:

 * `-h`   ~ Print usage information.
 * `-t`   ~ The time to be simulated, default value is 200, the unit is seconds.
 * `-c`   ~ The number of phases used. After every phases, the grid and the
            resources are adapted. The default value is 200, but the benchmark results
            presented in this thesis were obtained using 10 phases.

The benchmark results are written in CSV format to the current working directory.

The benchmarks have to be started with `mpirun`, the number of processes must be 
sufficient for the benchmarks. If the program terminates with an error from `libmpidynres`
right after starting it, try increasing the number of processes. The benchmarks are optimized
for using them with 56 processes. If there are not enough cores available, specify the 
`--oversubscribe` flag for `mpirun`.

For example, the benchmark for using `libmpidynres` on a Split by Four Grid
can be run with

`mpirun -n 56 --oversubscribe ./applications/build/SWE_p4est_benchMpidynresSplit_release -c 10`


### Examples ###

The example simulations accept the following parameters:

 * `-h, --help`                ~ Print usage information.
 * `-o, --output-basepath`     ~ Output base file name, mandatory argument, no default.
 * `-s, --scenario`            ~ Load scenario, (`artificial` | `dambreak`), default: `dambreak`
 * `-b, --bathymetry-path`     ~ Load bathymetry input from a file, overwrites `-s` option
 * `-d, --displacement-path`   ~ Load displacement data from a file, to be used in combination with `-b`
 * `-t, --simulated-time`      ~ Time in seconds to be simulated, default: 200
 * `-c, --num-phases`          ~ Number of phases, default: 200

The benchmarks have to be started with `mpirun`, the number of processes must be 
sufficient for the benchmarks. If the program terminates with an error from `libmpidynres`
right after starting it, try increasing the number of processes. If there are not enough 
cores available, specify the `--oversubscribe` flag for `mpirun`.

An example execution for the two dimensional simulation for 100s with 50 phases
would be started by calling

`mpirun -n 16 --oversubscribe ./applications/build/SWE_p4est_mpidynres2d_release -t 100 -c 50 -o output_file`
