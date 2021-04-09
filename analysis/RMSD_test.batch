#!/bin/bash	

#BSUB -P csc299
#BSUB -W 1:00
#BSUB -nnodes 1
#BSUB -J RMSD
#BSUB -o RMSD_%J.out
#BSUB -e RMSD_%J.err

hostname

date

which python
which spack
which gcc
which nvcc
which bpls
which cmake
which mpicc
which mpirun

cd /gpfs/alpine/scratch/iyakushin/csc299/Test3/entk_cvae_md_devel3/entk_cvae_md/analysis

jsrun -n 1 -c 42 -g 6 -a 1 -b packed:42 -d packed  python adios2pandasR4.py ../../../aggregate/aggregator0.bp 42 > 0.out 2>0.err


date
