#!/bin/bash

#SBATCH -J case2
#SBATCH -p second
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --ntasks-per-node=48
#SBATCH --mem=500G
# #SBATCH --exclude=tms09
#SBATCH --nodelist=tms12
#SBATCH -o d_case2.o%j
#SBATCH -e d_case2.e%j
# #SBATCH --exclusive

ulimit -s unlimited
# ulimit -n 4096
# export LD_LIBRARY_PATH=/opt/netcdf-c-4.9.2-fortran-4.6.0-gnu/lib:/opt/hdf5-1.14.0-gnu/lib:$LD_LIBRARY_PATH
# source /home/xiangjy23/.bashrc
mpirun -np 3 run/mksrfdata.x   run/case2.nml
wait

mpirun -np 3 run/mkinidata.x   run/case2.nml
wait

mpirun -np 3 run/colm.x       run/case2.nml
