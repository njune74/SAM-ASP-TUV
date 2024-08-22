#!/bin/sh
#$ -S /bin/csh
#$ -N W1101
#$ -o qsub_W1101.out
cd /rotor/data/p2033/SAM-ASP_ozone/
echo W1101 started at `date`
echo W1101 > /rotor/data/p2033/SAM-ASP_ozone/CaseName
/opt/sgi/mpt/mpt-2.02/bin/mpirun -np 12 ./SAM_RAD_CAM_MICRO_M2005_w_aero /rotor/data/p2033/SAM-ASP_ozone/RUNS/W1101/ > /rotor/data/p2033/SAM-ASP_ozone/OUT_RUNTIME/W1101.out 
echo W1101 ended at `date`
