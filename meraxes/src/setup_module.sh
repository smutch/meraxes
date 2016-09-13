module unload openmpi fftw hdf5 gsl
module load hdf5/x86_64/gnu/1.10.0-openmpi-1.10.2-psm
module load gsl/x86_64/gnu/1.9
module load openmpi/x86_64/intel/1.10.2-psm 
export PATH="/home/smutch/3rd_party/fftw-3.3.3/bin:$PATH"                       
export LIBRARY_PATH="/home/smutch/3rd_party/fftw-3.3.3/lib:$LIBRARY_PATH"       
export LD_LIBRARY_PATH="/home/smutch/3rd_party/fftw-3.3.3/lib:$LD_LIBRARY_PATH" 
export C_INCLUDE_PATH="/home/smutch/3rd_party/fftw-3.3.3/include:$C_INCLUDE_PATH"
