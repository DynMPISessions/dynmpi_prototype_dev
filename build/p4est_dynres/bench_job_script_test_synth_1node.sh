#!/bin/bash 


echo "running test nr. 1 with parameters: -np 4 -c 200 -l 32  -m s+ -n 4 -f 10 -b 1"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_b

echo "running test nr. 2 with parameters: -np 4 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4_nb

echo "running test nr. 3 with parameters: -np 8 -c 200 -l 32  -m s+ -n 4 -f 10 -b 1"
prterun -np 8 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_b

echo "running test nr. 4 with parameters: -np 8 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 8 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart8_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart8_nb

echo "running test nr. 5 with parameters: -np 12 -c 200 -l 32  -m s+ -n 4 -f 10 -b 1"
prterun -np 12 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_b

echo "running test nr. 6 with parameters: -np 12 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 12 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart12_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart12_nb

echo "running test nr. 7 with parameters: -np 16 -c 200 -l 32  -m s+ -n 4 -f 10 -b 1"
prterun -np 16 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16_b

echo "running test nr. 8 with parameters: -np 16 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 16 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart16_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart16_nb

echo "running test nr. 9 with parameters: -np 20 -c 200 -l 32  -m s+ -n 4 -f 10 -b 1"
prterun -np 20 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart20_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart20_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart20_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart20_b

echo "running test nr. 10 with parameters: -np 20 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 20 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart20_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart20_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart20_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart20_nb

echo "running test nr. 11 with parameters: -np 24 -c 200 -l 32  -m s+ -n 4 -f 10 -b 1"
prterun -np 24 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart24_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart24_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart24_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart24_b

echo "running test nr. 12 with parameters: -np 24 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 24 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart24_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart24_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart24_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart24_nb

echo "running test nr. 13 with parameters: -np 28 -c 200 -l 32  -m s+ -n 4 -f 10 -b 1"
prterun -np 28 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_b

echo "running test nr. 14 with parameters: -np 28 -c 200 -l 32  -m s+ -n 4 -f 10 -b 0"
prterun -np 28 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m s+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart28_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart28_nb

echo "running test nr. 15 with parameters: -np 4 -c 120 -l 32  -m i+ -n 4 -f 10 -b 1"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m i+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4inc_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_b

echo "running test nr. 16 with parameters: -np 4 -c 120 -l 32  -m i+ -n 4 -f 10 -b 0"
prterun -np 4 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 32  -m i+ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4inc_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4inc_nb

echo "running test nr. 17 with parameters: -np 32 -c 120 -l 4  -m i_ -n 4 -f 10 -b 1"
prterun -np 32 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 4  -m i_ -n 4 -f 10 -b 1
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4dec_b
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4dec_b
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4dec_b

echo "running test nr. 18 with parameters: -np 4 -c 120 -l 32  -m i+ -n 4 -f 10 -b 0"
prterun -np 32 --mca btl_tcp_if_include eth0 -H root-slurm-node01:40 -x LD_LIBRARY_PATH -x DYNMPI_BASE /opt/hpc/build/p4est_dynres/applications/build/SWE_p4est_benchOmpidynresSynthetic_debug -c 200 -l 4  -m i_ -n 4 -f 10 -b 0
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/synth/nstart4dec_nb
mkdir ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4dec_nb
mv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/*.csv ${DYNMPI_BASE}/build/p4est_dynres/applications/output/prrte/nstart4dec_nb


exit 0

