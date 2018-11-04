#!/bin/sh
#SBATCH -t 00:01:00 # execution time hh:mm:ss *OB*
#SBATCH -n 16 #tasks (MPI processes)
#SBATCH -c 24 #cores/task (shared-mem threads/process)
##SBATCH -N 16 #nodes
#SBATCH -p cola-corta
source env.sh
module load gcc/6.4.0
module load openmpi/2.1.1 extrae/3.5.2
export TRACE_NAME=trace_filename.prv # *OP*
srun -n 16 ${HOME}/trace.sh ./dgesv 2048
