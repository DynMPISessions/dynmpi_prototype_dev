#!/bin/bash 


echo "running test nr. 1 with parameters: -np 2 -c 200 -l 16  -m s+ -n 2 -f 10 -b 1"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2_b

echo "running test nr. 2 with parameters: -np 2 -c 200 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2_nb

echo "running test nr. 3 with parameters: -np 4 -c 200 -l 16  -m s+ -n 2 -f 10 -b 1"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_b

echo "running test nr. 4 with parameters: -np 4 -c 200 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_nb

echo "running test nr. 5 with parameters: -np 6 -c 200 -l 16  -m s+ -n 2 -f 10 -b 1"
prterun -np 6 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart6_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart6_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart6_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart6_b

echo "running test nr. 6 with parameters: -np 6 -c 200 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 6 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart6_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart6_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart6_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart6_nb

echo "running test nr. 7 with parameters: -np 8 -c 200 -l 16  -m s+ -n 2 -f 10 -b 1"
prterun -np 8 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_b

echo "running test nr. 8 with parameters: -np 8 -c 200 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 8 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_nb

echo "running test nr. 9 with parameters: -np 10 -c 200 -l 16  -m s+ -n 2 -f 10 -b 1"
prterun -np 10 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart10_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart10_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart10_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart10_b

echo "running test nr. 10 with parameters: -np 10 -c 200 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 10 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart10_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart10_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart10_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart10_nb

echo "running test nr. 11 with parameters: -np 12 -c 200 -l 16  -m s+ -n 2 -f 10 -b 1"
prterun -np 12 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_b

echo "running test nr. 12 with parameters: -np 12 -c 200 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 12 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_nb

echo "running test nr. 13 with parameters: -np 14 -c 200 -l 16  -m s+ -n 2 -f 10 -b 1"
prterun -np 14 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_b

echo "running test nr. 14 with parameters: -np 14 -c 200 -l 16  -m s+ -n 2 -f 10 -b 0"
prterun -np 14 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 200 -l 16  -m s+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_nb

echo "running test nr. 15 with parameters: -np 2 -c 120 -l 16  -m i+ -n 2 -f 10 -b 1"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 120 -l 16  -m i+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2inc_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2inc_b

echo "running test nr. 16 with parameters: -np 2 -c 120 -l 16  -m i+ -n 2 -f 10 -b 0"
prterun -np 2 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 120 -l 16  -m i+ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart2inc_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart2inc_nb

echo "running test nr. 17 with parameters: -np 16 -c 120 -l 2 -m i_ -n 2 -f 10 -b 1"
prterun -np 16 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 120 -l 2  -m i_ -n 2 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16dec_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16dec_b

echo "running test nr. 18 with parameters: -np 16 -c 120 -l 16  -m i+ -n 2 -f 10 -b 0"
prterun -np 16 --mca btl_tcp_if_include eth0 -H root-node01:2,root-node02:2,root-node03:2,root-node04:2,root-node05:2,root-node06:2,root-node07:2,root-node08:2, -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_release -c 120 -l 2  -m i_ -n 2 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16dec_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16dec_nb


exit 0

