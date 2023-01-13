#!/bin/bash

#SBATCH --job-name=diablo3_test
#
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=15G
#SBATCH --partition=serc
#SBATCH --constraint=[CLASS:SH3_CBASE|CLASS:SH3_CPERF]
#SBATCH --output=output.log
#SBATCH --error=error.log
#SBATCH --mail-user=hilditch@stanford.edu
#SBATCH --mail-type=end,fail

ml load hdf5/1.12.2
ml load fftw/3.3.10
# ml load valgrind

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/users/hilditch/packages/toml-f/lib64
# export OMPI_MCA_pml=^ucx


$GROUP_HOME/hilditch/scripts/diablo3_test/for/setup_run 8 2 --toml > setup.log

srun -n 16 ./diablo > diablo.log
# srun -n 32 valgrind --log-file=valgrind.log --track-origins=yes --error-limit=no --leak-check=full --show-leak-kinds=all --suppressions=/share/software/user/open/openmpi/4.1.2/share/openmpi/openmpi-valgrind.supp ./diablo > diablo.log

