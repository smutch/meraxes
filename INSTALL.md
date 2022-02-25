**THIS IS A WORK IN PROGRESS!**

# Dependencies

The following instructions list possible methods to get and install the required dependencies for Meraxes. These are intended for development use only. If actually running on an HPC facility then all required libraries should already be available to you. Please consult your system's help pages or contact your sys-admin for further information.

## All platforms

A simple (although heavy) method of aquiring all the necessary dependencies is to use [Spack](https://spack.readthedocs.io/en/latest/). A `spack.yml` config file is provided in the Meraxes repo. To use this method, first ensure you have spack installed and activated in your current shell instance. Then simply run the provided `bootstrap-dev.sh` script and follow the prompts.

## Linux / Unix

### CMAKE, GSL, HDF5

If you have root access, then installing via your distribution's package manager is your best bet.  If you are working on a supercomputer then chances these are all already installed or are available as modules.  Please note that HDF5 must be compiled with MPI bindings.  If you're having problems then please contact your local sysadmin.

### FFTW3

Although fftw3 is likely also available, we require the `--enable-mpi` and `--enable-float` flags which are somewhat non-standard.  If you need to download and install fftw the old fashioned way then following something like the below steps should do the trick.  Please note that you should replace the exact version number with whatever is the most current one listed on the [fftw website](http://www.fftw.org/).

``` sh
mkdir -p ~/3rd_party/fftw3 && cd !$
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar xzf fftw-3.3.10.tar.gz
mv fftw-3.3.10 src && cd !$
CC=`which mpicc` ./configure --prefix=$HOME/3rd_party/fftw3 --enable-mpi --enable-float
make && make install && make clean
```

## Mac

By far the easiest way to install 3rd party libraries on a mac is via [homebrew](https://brew.sh)...

### CMAKE, GSL & HDF5

```sh
brew install cmake
brew install gsl
brew install hdf5 --with-mpi
```

### FFTW3

Although fftw3 can be installed via homebrew, we require the `--enable-mpi` and `--enable-float` flags which are currently not an option with `brew install fftw`.  Instead, we will need to do it the good old fashioned way.

``` sh
mkdir -p ~/3rd_party/fftw3 && cd !$
wget http://www.fftw.org/fftw-3.3.10.tar.gz
tar xzf fftw-3.3.10.tar.gz
mv fftw-3.3.10 src && cd !$
CC=`which mpicc` ./configure --prefix=$HOME/3rd_party/fftw3 --enable-mpi --enable-float
make && make install && make clean
```

# Compiling Meraxes

There are a few ways to do this, but the only way which is independent of the cmake version being used is the following:

```sh
cd src
mkdir build
cd build
cmake ..
make
```
