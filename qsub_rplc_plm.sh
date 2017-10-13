#!/bin/bash --login

# mpirun -np 6 lmp_rplc_plm -p 6x1 -in in.replica -screen none
mpirun -np 2 lmp_rplc_plm -p 2x1 -in in.rplc_plm 
