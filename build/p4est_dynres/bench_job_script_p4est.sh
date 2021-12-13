#!/bin/bash 

#SBATCH --job-name=bench_job
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=domi-huber@tum.de
#SBATCH --ntasks=28
#SBATCH --mem=1gb
#SBATCH --time=00:50:00
#SBATCH --output=bench_job_output.log

pwd; hostname; date



export LD_LIBRARY_PATH=${DYNMPI_BASE}/build/p4est_dynres/p4est/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${DYNMPI_BASE}/build/p4est_dynres/libmpidynres/build/lib:$LD_LIBRARY_PATH

#build benchmark
cd ${DYNMPI_BASE}/build/p4est_dynres/applications
scons examples=benchOmpidynresSynthetic compileMode=release

#Generate host file
cd ${DYNMPI_BASE}
scontrol show hostname ${SLURM_JOB_NODELIST} > bench_base_folder/hostfile-${SLURM_JOB_ID}.txt

echo "P4EST: running test nr. 1 with parameters: -np 28 -c 16 -l 224  -m i+ -n 2 -f 1 -b 1"
prterun -np 28 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 224 -t 200 -m i+ -n 2 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_b


echo "P4EST: running test nr. 2 with parameters: -np 28 -c 16 -l 224  -m i+ -n 2 -f 10 -b 0"
prterun -np 28 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 224 -t 200 -m i+ -n 2 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_nb


echo "P4EST: running test nr. 3 with parameters: -np 224 -c 16 -l 28 -m i_ -n 2 -f 1 -b 1"
prterun -np 224 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 28 -t 200 -m i_ -n 2 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart16dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2dec_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart16dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2dec_b


echo "P4EST: running test nr. 4 with parameters: -np 224 -c 16 -l 28 -m i_ -n 2 -f 10 -b 0"
prterun -np 224 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 28 -t 200 -m i_ -n 2 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart16dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2dec_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart16dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2dec_nb


echo "P4EST: running test nr. 5 with parameters: -np 28 -c 16 -l 224  -m s+ -n 2 -f 1 -b 1"
prterun -np 28 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 224 -t 200 -m s+ -n 2 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_b


echo "P4EST: running test nr. 6 with parameters: -np 28 -c 20 -l 224  -m s+ -n 2 -f 10 -b 0"
prterun -np 28 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 24 -l 224 -t 200 -m s+ -n 2 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_nb


exit 0

