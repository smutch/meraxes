**THIS IS A WORK IN PROGRESS!**

# Dependencies

## Linux / Unix

### GSL, HDF5 & Scons

If you have root access, then installing via your distribution's package manager is your best bet.  If you are working on a supercomputer then chances are GSL and HDF5 are already installed or are available as modules.  Please note that HDF5 must be compiled with MPI bindings.  If you're having problems then please contact your local sysadmin.

### FFTW3

Although fftw3 is likely also available, we require the `--enable-mpi` and `--enable-float` flags which are often not standard.  If you need to download and install fftw the old fashioned way then following something like the below steps should do the trick.  Please note that you should replace the exact version number with whatever is the most current one listed on the [fftw website](http://www.fftw.org/).

``` sh
mkdir -p ~/3rd_party/fftw3 && cd !$
wget http://www.fftw.org/fftw-3.3.6-pl2.tar.gz
tar xzf fftw-3.3.6-pl2.tar.gz
mv fftw-3.3.6-pl2 src && cd !$
CC=`which mpicc` ./configure --prefix=$HOME/3rd_party/fftw3 --enable-mpi --enable-float
make && make install && make clean
```

## Mac

By far the easiest way to install 3rd party libraries on a mac is via [homebrew](https://brew.sh)...

### GSL, HDF5 & Scons

```sh
brew install gsl
brew install hdf5 --with-mpi
brew install scons
```

### FFTW3

Although fftw3 can be installed via homebrew, we require the `--enable-mpi` and `--enable-float` flags which are currently not an option with `brew install fftw`.  Instead, we will need to do it the good old fashioned way.

``` sh
mkdir -p ~/3rd_party/fftw3 && cd !$
wget http://www.fftw.org/fftw-3.3.6-pl2.tar.gz
tar xzf fftw-3.3.6-pl2.tar.gz
mv fftw-3.3.6-pl2 src && cd !$
CC=`which mpicc` ./configure --prefix=$HOME/3rd_party/fftw3 --enable-mpi --enable-float
make && make install && make clean
```

# Compiling Meraxes

In order to compile Meraxes, you need to first set the paths to all of the dependencies.  This is done using the files in the `src/deps` directory.  Either modify the paths in the file corresponding to your system, or create a new file in the `src/deps` directory called `custom.py` and populate it appropriately. 

After that, the code can be compiled by running `scons` in the `src` directory.
