#!/bin/sh
#SBATCH --job-name=SND1E1S1_v10
#SBATCH --nodes 1
#SBATCH --ntasks 16 
#SBATCH --mail-user=nicole.june@colostate.edu
#SBATCH --mail-type=END

cd /pierce-scratch/njune/SAM-ASP_tuv/
echo SND1 started at `date`
echo SND1 > /pierce-scratch/njune/SAM-ASP_tuv/CaseName

##from rotor## LD_LIBRARY_PATH=/rotor/data/p1808/CMAQv5.0.1/lib/x86_64/ifort/mpich/lib:${LD_LIBRARY_PATH}
LD_LIBRARY_PATH=/usr/local/mpich/3.2-intel18.0.0/lib:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH >> /pierce-scratch/njune/SAM-ASP_tuv/log.txt

#AD DEBUG. Looks like the old case was run with /opt/sgi/mpt/mpt-2.02/bin/mpirun -np 12
#try that instead of mpich?
mpirun -np 16 ./SAM_RAD_CAM_MICRO_M2005_w_aero_SND1E1S1_v10 /pierce-scratch/njune/SAM-ASP_tuv/RUNS/SND1/ > /pierce-scratch/njune/SAM-ASP_tuv/OUT_RUNTIME/SND1E1S1_v10.out   
  
echo SND1 ended at `date`
#/opt/sgi/mpt/mpt-2.02/bin/mpirun
