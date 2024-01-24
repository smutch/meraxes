# BUILDING MERAXES

Meraxes has been compiled and run successfully on the following operating systems and architectures:

- Mac OS X >=10.10 (Yosemite) on x86_64
- Mac OS X >=10.13 (Sierra) on aarch64
- Linux (Ubuntu 16.04, 18.04) on x86_64
- Linux (CentOS 7) on x86_64

We expect things to work without any problem on any modern Linux or Mac system.  If you have any problems then please open an issue.

## Dependencies

Meraxes has a few dependencies which must be installed before you can compile it.  These are:

- [CMAKE](https://cmake.org)
- [GSL](https://www.gnu.org/software/gsl/)
- [HDF5](https://www.hdfgroup.org)
- [FFTW3](http://www.fftw.org)

For Nix/NixOS users, a flake is provided to install all dependencies in a dev shell.
For Spack users, a dev environment is provided in `spack-env`.
For everyone else, please see the system-dependent instructions below.

### Linux / Unix

#### CMAKE, GSL, HDF5

If you have root access, then installing via your distribution's package manager is your best bet.  If you are working on a supercomputer then chances these are all already installed or are available as modules.  Please note that HDF5 must be compiled with MPI bindings.  If you're having problems then please contact your local sysadmin.

#### FFTW3

Although fftw3 is likely also available, we require the `--enable-mpi` and `--enable-float` flags which are somewhat non-standard.  If you need to download and install fftw the old fashioned way then following something like the below steps should do the trick.  Please note that you should replace the exact version number with whatever is the most current one listed on the [fftw website](http://www.fftw.org/).

``` sh
mkdir -p ~/3rd_party/fftw3 && cd !$
wget http://www.fftw.org/fftw-3.3.6-pl2.tar.gz
tar xzf fftw-3.3.6-pl2.tar.gz
mv fftw-3.3.6-pl2 src && cd !$
CC=`which mpicc` ./configure --prefix=$HOME/3rd_party/fftw3 --enable-mpi --enable-float
make && make install && make clean
```

### Mac

By far the easiest way to install 3rd party libraries on a mac is via [homebrew](https://brew.sh). The only problem here is that we have no way to pin versions of the libraries we need so guaranteeing things you compile now will work in the future is not possible with our current build set up. You might want to consider Nix or Spack for more a more dependable dev environment.

#### CMAKE, GSL & HDF5

```sh
brew install cmake
brew install gsl
brew install hdf5 --with-mpi
```

#### FFTW3

Although fftw3 can be installed via homebrew, we require the `--enable-mpi` and `--enable-float` flags which are currently not an option with `brew install fftw`.  Instead, we will need to do it the good old fashioned way.

``` sh
mkdir -p ~/3rd_party/fftw3 && cd !$
wget http://www.fftw.org/fftw-3.3.6-pl2.tar.gz
tar xzf fftw-3.3.6-pl2.tar.gz
mv fftw-3.3.6-pl2 src && cd !$
CC=`which mpicc` ./configure --prefix=$HOME/3rd_party/fftw3 --enable-mpi --enable-float
make && make install && make clean
```

## Compiling

There are a few ways to do this, but the only way which is independent of the cmake version being used is the following:

```sh
mkdir build
cd build
cmake ..
make
```

Meraxes has a number of compile-time options and parameters:

BUILD_SHARED_LIBS
: Build Meraxes as a shared library. Default is OFF.

BUILD_TESTS
: Build the test suite. Default is OFF.

CALC_MAGS
: Calculate magnitudes. Default is OFF.

SECTOR_ROOT
: Path to the sector repo. Default is `src/sector`. This is only used if `CALC_MAGS=ON`.

MAGS_N_BANDS
: Number of bands to calculate magnitudes for. Only applicable if `CALC_MAGS=ON`. Default is 6.

MAGS_N_SNAPS
: Number of snapshots to calculate magnitudes for. Only applicable if `CALC_MAGS=ON`. Default is 3.

CMAKE_BUILD_TYPE
: Build type. Default is RelWithDebInfo (Release with debug symbols).

CMAKE_INSTALL_PREFIX
: Install prefix. Default is a directory called `target` in the repo base directory.

ENABLE_PROFILING
: Enable profiling. Default is OFF.

GDB
: Enable GDB debugging. Default is OFF.

N_HISTORY_SNAPS
: Number of snapshots to star formation history for. This is simulation dependent. If the value is too low then Meraxes will crash with an error telling you what the minimum value is for your simulation. Default is 5.

USE_CUDA
: Use CUDA accelerated reionization calculation. Default is OFF.

USE_MINI_HALOS
: Use mini-halos. Default is ON.

You can set these on the command line when running cmake, e.g.:

```sh
cmake .. -DN_HISTORY_SNAPS=10 -DCMAKE_INSTALL_PREFIX=/path/to/install
```

Alternatively, you can use the ccmake TUI to set these options.  This is useful if you want to see what the default values are.

```sh
ccmake ..
```

Lastly, you can specify any options in a `local.cmake` file.  This is useful if you want to set some options permanently.  For example, if you want to always build Meraxes as a shared library then you can create a `local.cmake` file with the following contents:

```cmake
set(BUILD_SHARED_LIBS ON)
```
