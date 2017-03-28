**THIS IS A WORK IN PROGRESS!**

# Mac

By far the easiest way to install 3rd party libraries on a mac is via [homebrew](https://brew.sh)...

## GSL

```sh
brew install gsl
```

## HDF5

``` sh
brew install hdf5 --with-mpi
```

## FFTW3

Although fftw3 can be installed via homebrew, we require the `--enable-mpi` and `--enable-float` flags which are currently not an option with `brew install fftw`.  Instead, we will need to do it the good old fashioned way.

``` sh
mkdir -p ~/3rd_party/fftw3 && cd !$
wget http://www.fftw.org/fftw-3.3.6-pl2.tar.gz
tar xzf fftw-3.3.6-pl2.tar.gz
mv fftw-3.3.6-pl2 src && cd !$
CC=`which mpicc` ./configure --prefix=$HOME/3rd_party/fftw3 --enable-mpi --enable-float
make && make install && make clean
```
