#!/bin/bash 

echo "P4EST: running test nr. 1 with parameters: -np 4 -c 200 -l 32  -m i+ -n 4 -f 1 -b 1"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_debug -c 16 -l 32 -t 1 -m i+ -n 4 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4inc_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_b

rm -R /dev/shm/*

echo "P4EST: running test nr. 2 with parameters: -np 4 -c 200 -l 32  -m i+ -n 4 -f 10 -b 0"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_debug -c 16 -l 32 -t 1 -m i+ -n 4 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4inc_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_nb

rm -R /dev/shm/*

echo "P4EST: running test nr. 3 with parameters: -np 32 -c 200 -l 4 -m i_ -n 4 -f 1 -b 1"
prterun -np 32 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_debug -c 16 -l 4 -t 1 -m i_ -n 4 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart32dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4dec_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart32dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4dec_b

rm -R /dev/shm/*

echo "P4EST: running test nr. 4 with parameters: -np 32 -c 200 -l 4 -m i_ -n 4 -f 10 -b 0"
prterun -np 32 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_debug -c 16 -l 4 -t 1 -m i_ -n 4 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart32dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4dec_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart32dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4dec_nb

rm -R /dev/shm/*

echo "P4EST: running test nr. 5 with parameters: -np 4 -c 200 -l 32  -m s+ -n 4 -f 1 -b 1"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_debug -c 16 -l 32 -t 1 -m s+ -n 4 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4seq_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4seq_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4seq_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4seq_b

rm -R /dev/shm/*

echo "P4EST: running test nr. 6 with parameters: -np 4 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_debug -c 16 -l 32 -t 1 -m s+ -n 4 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4seq_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart4seq_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4seq_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4seq_nb

rm -R /dev/shm/*

exit 0

