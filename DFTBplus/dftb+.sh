#! /bin/bash

module purge 
module load intel impi
module load dftb+

export OMP_NUM_THREADS=$SLURM_NPROCS

python run_dftb+.py

srun dftb+

exit 0
