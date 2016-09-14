#!/bin/bash
module purge
module load openmpi/x86_64/intel/1.10.2-psm gsl/x86_64/gnu/1.9  hdf5/x86_64/gnu/1.10.0-openmpi-1.10.2-psm 
#for libp in `grep 'libp' deps/gstar.py| awk 'BEGIN{FS=":"}{print $2}'`; do export LD_LIBRARY_PATH=`echo $libp| tr -d "'"|tr -d ','`:$LD_LIBRARY_PATH; done
