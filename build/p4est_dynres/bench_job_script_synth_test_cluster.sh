#!/bin/bash 

#SBATCH --job-name=bench_job
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=domi-huber@tum.de
#SBATCH --ntasks=56
#SBATCH --mem=1gb
#SBATCH --time=00:20:00
#SBATCH --output=bench_job_output.log

pwd; hostname; date

#Generate host file
cd ${DYNMPI_BASE}
scontrol show hostname ${SLURM_JOB_NODELIST} > ${DYNMPI_BASE}/hostfile-${SLURM_JOB_ID}.txt


echo "running test nr. 1 with parameters: -np 28 -c 120 -l 56 -m i+ -n 28 -f 10 -b 1"
prterun -np 28 --mca btl_tcp_if_include eth0 --hostfile ${DYNMPI_BASE}/hostfile-${SLURM_JOB_ID}.txt -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 120 -l 56  -m i+ -n 28 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2inc_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2inc_b


exit 0

