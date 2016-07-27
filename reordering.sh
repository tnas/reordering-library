#$ -S /bin/bash
#$ -q fila
#$ -pe threads 8

export OMP_NUM_THREADS=8

./reordering-library -t 1

# Script for executing the reordering-library program into 
# the cluster Altix-XE from LNCC.
# Commando for submit the program to the cluster's queue: 
# qsub reordering.sh