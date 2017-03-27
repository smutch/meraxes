# Mac library dependencies
import os

deps = {}

deps['exec'] = {
    'mpicc' : 'mpicc',
    'cc' : 'cc',
    'git' : 'git',
}

deps['gsl'] = {
    'inclp' : '/usr/local/Cellar/gsl/2.3/include',
    'libp'  : '/usr/local/Cellar/gsl/2.3/lib',
    'lib'   : ['gsl', 'gslcblas'],
}

deps['hdf5'] = {
    'inclp' : '/usr/local/Cellar/hdf5/1.8.16_1/include',
    'libp'  : '/usr/local/Cellar/hdf5/1.8.16_1/lib',
    'lib'   : ['hdf5', 'hdf5_hl', 'z'],
}

deps['gbpCode'] = {
    'inclp' : os.path.expanduser('~/3rd_party/gbpCode/myInclude'),
    'libp'  : os.path.expanduser('~/3rd_party/gbpCode/myLib/mpi'),
    'lib'   : 'gbpLib',
}

deps['fftw'] = {
    'inclp' : os.path.expanduser('~/3rd_party/fftw-3.3.5/include'),
    'libp'  : os.path.expanduser('~/3rd_party/fftw-3.3.5/lib'),
    'lib'   : ['fftw3f_mpi', 'fftw3f'],
}

deps['mlog'] = {
    'inclp' : os.path.join(os.path.dirname(os.path.realpath(__file__)), '../mlog'),
    'libp' : os.path.join(os.path.dirname(os.path.realpath(__file__)), '../mlog'),
    'lib'   : ['mlog'],
}
