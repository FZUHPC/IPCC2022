#!/bin/bash
#SBATCH -p amd_256
#SBATCH -J pivot
#SBATCH -t 1
#SBATCH -N 1
#SBATCH -o result.log
#SBATCH -e error.log
source /public1/soft/modules/module.sh 
export OMP_PROC_BIND=TRUE
g++ -std=c++11 -O3  -fopenmp -mfma  -ftree-vectorize -ffast-math  -funroll-all-loops -mavx2 -mtune=native -march=native pivot.c -o pivot 
srun -p amd_256 -N 1 ./pivot

