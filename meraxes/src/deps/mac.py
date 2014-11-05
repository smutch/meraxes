# Mac library dependencies
import os

deps = {}

deps['exec'] = {
    'mpicc' : 'mpicc',
    'cc' : 'cc',
    'git' : 'git',
}

deps['gsl'] = {
    'inclp' : None,
    'libp'  : None,
    'lib'   : ['gsl', 'gslcblas'],
}

deps['hdf5'] = {
    'inclp' : None,
    'libp'  : None,
    'lib'   : ['hdf5', 'hdf5_hl', 'z'],
}

deps['gbpCode'] = {
    'inclp' : os.path.expanduser('~/3rd_party/gbpCode/myInclude'),
    'libp'  : os.path.expanduser('~/3rd_party/gbpCode/myLib'),
    'lib'   : 'gbpLib',
}

deps['fftw'] = {
    'inclp' : os.path.expanduser('~/3rd_party/fftw-3.3.3/include'),
    'libp'  : os.path.expanduser('~/3rd_party/fftw-3.3.3/lib'),
    'lib'   : ['fftw3f_omp', 'fftw3f', 'gomp'],
}

deps['21cmfast'] = {
    'inclp' :'../../21cmfast/include',
    'libp'  :'../../21cmfast/lib',
    'lib'   : ['21cmfast', 'gomp'],
}

