#!/bin/sh
#$ -S /bin/csh
#$ -N tuvtest
#$ -o qsub_tuvtest.out
cd /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/
echo tuvtest started at `date` >> /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/log.txt
echo tuvtest > /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/CaseName

##from rotor## LD_LIBRARY_PATH=/rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/mpich/lib:${LD_LIBRARY_PATH}
LD_LIBRARY_PATH=/usr/lib64/mpich-3.2/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH >> /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/log.txt

#AD DEBUG. Looks like the old case was run with /opt/sgi/mpt/mpt-2.02/bin/mpirun -np 12
#try that instead of mpich?
mpirun -np 1 ./SAM_RAD_CAM_MICRO_M2005_w_aero /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/RUNS/tuvtest/ >& /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/OUT_RUNTIME/tuvtest-test33-MM.out  
  
echo tuvtest ended at `date` >> /nas/data_cellar/rotor/data/p2033/SAM-ASP_tuv/log.txt
#/opt/sgi/mpt/mpt-2.02/bin/mpirun
