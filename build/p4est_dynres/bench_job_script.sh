#!/bin/bash 

#SBATCH --job-name=bench_job
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=domi-huber@tum.de
#SBATCH --ntasks=28
#SBATCH --mem=1gb
#SBATCH --time=00:50:00
#SBATCH --output=bench_job_output.log

pwd; hostname; date

#load required modules
module unload intel-mpi
module load zlib
module load libtool
module load hwloc
module load m4
module load autoconf
module load automake
module load flex
module load hdf5
module load netcfd

#build own mpi,pmix,prrte,p4est,libmpidynres
cd ${DYNMPI_BASE}/build/openpmix
./../bin/build-openpmix.sh
cd ${DYNMPI_BASE}/build/prrte
./../bin/build-prrte.sh
cd ${DYNMPI_BASE}/build/ompi
./../bin/build-openmpi.sh
cd ${DYNMPI_BASE}/build/p4est_dynres/p4est
./configure --enable-mpi --without-blas && make && make install
cd ${DYNMPI_BASE}/build/p4est_dynres/libmpidynres
make && make install

export LD_LIBRARY_PATH=${DYNMPI_BASE}/build/p4est_dynres/p4est/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${DYNMPI_BASE}/build/p4est_dynres/libmpidynres/build/lib:$LD_LIBRARY_PATH

#build benchmark
cd ${DYNMPI_BASE}/build/p4est_dynres/applications
scons examples=benchOmpidynresSynthetic compileMode=release

#Generate host file
cd ${DYNMPI_BASE}
scontrol show hostname ${SLURM_JOB_NODELIST} > bench_base_folder/hostfile-${SLURM_JOB_ID}.txt

echo "running test nr. 1 with parameters: -np 2 -c 400 -l 28  -m s+ -n 4 -f 40"

prterun -np 4 --mca btl_tcp_if_include eth0 -hostfile ${DYNMPI_BASE}/hostfile-${SLURM_JOB_ID}.txt -x LD_LIBRARY_PATH -x DYNMPI_BASE .${DYNMPI_BASE}/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 400 -l 28  -m i+ -n 4 -f 40

echo "test nr. 1 finished"

exit 0

