#!/bin/sh
source /opt/cesga/easybuild-cesga/software/MPI/gcc/6.4.0/openmpi/2.1.1/extrae/3.5.2/etc/extrae.sh
export EXTRAE_CONFIG_FILE=${HOME}/extrae1.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so # C code
$* # Run the program
