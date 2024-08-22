#!/bin/sh
#$ -S /bin/csh
#$ -N tuvtest
#$ -o qsub_tuvtest.out
cd /rotor/data/p2033/SAM-ASP_tuv/
echo tuvtest started at `date` >> /rotor/data/p2033/SAM-ASP_tuv/log.txt
echo tuvtest > /rotor/data/p2033/SAM-ASP_tuv/CaseName

LD_LIBRARY_PATH=/rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/mpich/lib:${LD_LIBRARY_PATH}
#LD_LIBRARY_PATH=/opt/sgi/mpt/mpt-2.02/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH >> /rotor/data/p2033/SAM-ASP_tuv/log.txt

#AD DEBUG. Looks like the old case was run with /opt/sgi/mpt/mpt-2.02/bin/mpirun -np 12
#try that instead of mpich?
/rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/mpich/bin/mpirun -np 1 ./SAM_RAD_CAM_MICRO_M2005_w_aero /rotor/data/p2033/SAM-ASP_tuv/RUNS/tuvtest/ >& /rotor/data/p2033/SAM-ASP_tuv/OUT_RUNTIME/tuvtest-test3-MM.out #AD DEBUG 


#/opt/sgi/mpt/mpt-2.02/bin/mpirun -np 12 ./SAM_RAD_CAM_MICRO_M2005_w_aero /rotor/data/p2033/SAM-ASP_tuv/RUNS/tuvtest/ >& /rotor/data/p2033/SAM-ASP_tuv/OUT_RUNTIME/tuvtest.out
  
echo tuvtest ended at `date` >> /rotor/data/p2033/SAM-ASP_tuv/log.txt
#/opt/sgi/mpt/mpt-2.02/bin/mpirun
