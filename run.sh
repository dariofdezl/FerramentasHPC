#!/bin/bash
#SBATCH -N 1 #(1 node)
#SBATCH -n 1 #(1 job)
#SBATCH -c 24 #(24 cores per job)
#SBATCH -t 10:00:00

source env.sh

make clean

make

make run

