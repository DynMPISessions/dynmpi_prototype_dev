# Towards Dynamic Resource Management with MPI Sessions and PMIx - Prototype Code for Benchmarks
## Introduction
This repository provides the code of the prototype developed in the context of my master's thesis "Towards Dynamic Resource Management with MPI Sessions and PMIx". 

The goal of this project is to develop resource adaptive MPI applications and runtime environments by extending the MPI Sessions model and PMIx. This prototype implementation is based on the Open-MPI code, i.e. the sessions_pr branch invloved in the following pull request: https://github.com/opne-mpi/ompi/pull/9097  


The prototype is currently under further development to extend its functionalities. 

## Prerequisits
Open-MPI:
* m4 1.4.17
* autoconf 2.69
* automake 1.15
* libtool 2.4.6
* flex 2.5.35
* hwloc 2.5
* libevent 2.1.12
* zlib (recommended)

Building the benchmarks:
* scons 

## Compiling and Installing
First export the path to the base directory (the directory including this README):

`export DYNMPI_BASE=[/path/to/base_dir]`

Export the ompi, open-pmix and prrte install paths:

`export PMIX_ROOT=$DYNMPI_BASE/install/pmix`

`export PRRTE_ROOT=$DYNMPI_BASE/install/prrte`

`export OMPI_ROOT=$DYNMPI_BASE/install/ompi`

Export the hwloc and libevent install paths:

`export HWLOC_INSTALL_PATH=[/path/to/hwloc]`

`export LIBEVENT_INSTALL_PATH=[/path/to/libevent]`

Building Open-PMIx (you might also want to consider the README.md in the build/openpmix directory):

`cd build/openpmix`

`./autogen.pl`

`mkdir build/openpmix/build`

`./../configure --prefix=$PMIX_ROOT -with-hwloc=$HWLOC_INSTALL_PATH --with-libevent=$LIBEVENT_INSTALL_PATH`

`make && make all install`

`cd ../..`

Building PRRTE (you might also need to consider the README.md in the build/prrte directory):

`cd build/prrte`

`./autogen.pl`

`./configure --prefix=$PRRTE_ROOT -with-hwloc=$HWLOC_INSTALL_PATH --with-libevent=$LIBEVENT_INSTALL_PATH -with-pmix=$PMIX_ROOT`

`make && make all install`

`cd ..`

Building Open-MPI (you might also need to consider the README.md in the build/ompi directory):

`cd build/ompi`

`./autogen.pl`

`./configure --prefix=$OMPI_ROOT -with-hwloc=$HWLOC_INSTALL_PATH --with-libevent=$LIBEVENT_INSTALL_PATH --with-pmix=$PMIX_ROOT --with-prrte=$PRRTE_ROOT --with-ucx=no --disable-mpi-fortran`

`make && make all install`

`cd ..`

Building P4est (you might also need to consider the README.md in the build/p4est_dynres/p4est directory):

`cd build/p4est_dynres/p4est`

(When compiling for the first time, it might be necessary to call `./bootstrap`)

`./configure --enable-mpi --without-blas`

`make && make install`

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DYNMPI_BASE/build/p4est_dynres/p4est/local/lib`

`cd ../..`

(Optional) Building libmpidynres (you might also need to consider the README.md in the build/p4est_dynres/libmpidynres directory):

`cd build/p4est_dynres/libmpidynres`

`make && make install`

`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DYNMPI_BASE/build/p4est_dynres/libmpidynres/build/lib`

`cd ../..`

Building the benchmarks:

`scons example=[benchOmpidynresSynthetic/benchOmpidynresFixed] compileMode=[debug/release]`

## Running the benchmarks
Note: The `DYNMPI_BASE environment variable has to be added to the environment on every invloved node. Using the -x option for the prterun command is mandatory but eventually not sufficient.

### Synthetic Benchmark:
* Incremental mode (Addition)

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_release -c 120 -l 122 -m i+ -n 28 -f 10 -b 0`

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_nb_release -c 120 -l 122 -m i+ -n 28 -f 10 -b 1`

* Incremental mode (Subtraction)

`prterun -np 122 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_release -c 120 -l 28 -m i_ -n 28 -f 10 -b 0`

`prterun -np 122 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_nb_release -c 120 -l 28 -m i_ -n 28 -f 10 -b 1`


* Step mode (Addition):

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_release -c 200 -l 122 -m s+ -n 28 -f 10 -b 0`

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_nb_release -c 200 -l 122 -m s+ -n 28 -f 10 -b 1`

* Step mode (Subtraction):

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_release -c 200 -l 122 -m s_ -n 28 -f 10 -b 0`

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_v1_nb_release -c 200 -l 122 -m s_ -n 28 -f 10 -b 1`

### SWE Benchmark:

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 10 -t 10 -l 122 -m i+ -n 28 -f 10 -b 1`

`prterun -np 28 --mca btl_tcp_if_include eth0 -H node01:28,node02:28,node03:28,node04:28 -x LD_LIBRARY_PATH -x DYNMPI_BASE $DYNMPI_BASE/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 10 -t 10 -l 122 -m i+ -n 28 -f 1 -b 0`

## Contact
In case of any questions please contact me via email: `domi.huber@tum.de`


