#!/usr/bin/bash -l
#SBATCH --time 1-0:0:0

module purge
module -q restore akantu

export OMP_NUM_THREADS=1
srun ./cohesive_extrinsic --aka_input_file ./material-elastic.dat
