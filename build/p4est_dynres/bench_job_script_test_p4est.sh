#!/bin/bash 

echo "P4EST: running test nr. 1 with parameters: -np 2 -c 16 -l 16  -m i+ -n 2 -f 1 -b 1"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 16 -t 1 -m i+ -n 2 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_b


echo "P4EST: running test nr. 2 with parameters: -np 2 -c 16 -l 16  -m i+ -n 2 -f 10 -b 0"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 16 -t 1 -m i+ -n 2 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2inc_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2inc_nb


echo "P4EST: running test nr. 3 with parameters: -np 16 -c 16 -l 2 -m i_ -n 2 -f 1 -b 1"
prterun -np 16 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 2 -t 1 -m i_ -n 2 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart16dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2dec_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart16dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2dec_b


echo "P4EST: running test nr. 4 with parameters: -np 16 -c 16 -l 2 -m i_ -n 2 -f 10 -b 0"
prterun -np 16 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 2 -t 1 -m i_ -n 2 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart16dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2dec_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart16dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2dec_nb


echo "P4EST: running test nr. 5 with parameters: -np 2 -c 16 -l 16  -m s+ -n 2 -f 1 -b 1"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 16 -l 16 -t 1 -m s+ -n 2 -f 1 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_b


echo "P4EST: running test nr. 6 with parameters: -np 2 -c 20 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresFixed_release -c 20 -l 16 -t 1 -m s+ -n 2 -f 1 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/p4est/nstart2seq_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/p4est_nstart2seq_nb


exit 0

