#!/bin/sh

cd /mnt/home/smutch/models/21cm_sam/meraxes/output/results/diagnostics
python ../allresults.py /lustre/projects/p070_astro/smutch/meraxes/tuning/output/meraxes.hdf5 --ext=png
python cold_gas_fraction.py /lustre/projects/p070_astro/smutch/meraxes/tuning/output/meraxes.hdf5 5
python gas_sd.py /lustre/projects/p070_astro/smutch/meraxes/tuning/output/meraxes.hdf5 5

