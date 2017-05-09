#!/bin/bash --login

#$ -q ibnet
#$ -pe orte 48
#$ -cwd 
#$ -j N
#$ -S /bin/bash
#$ -N finereplica

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/newhome/jhaberstroh/Downloads/matheval-compile/BUILD/lib/

mpirun -np ${NSLOTS} lmp_replica -p 6x8 -in in.replica -screen none
