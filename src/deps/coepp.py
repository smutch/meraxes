# COEPP library dependencies
import os

deps = {}

deps['exec'] = {
    'mpicc' : '/home/borg/smutch/cephfs/3rd_party/openmpi-2.0.1/bin/mpicc',
    'cc' : '/usr/bin/gcc',
    'git' : '/usr/bin/git',
}

gsl_dir = '/home/borg/smutch/cephfs/3rd_party/gsl-2.2.1'
deps['gsl'] = {
    'inclp' : os.path.join(gsl_dir, 'include'),
    'libp'  : os.path.join(gsl_dir, 'lib'),
    'lib'   : ['gsl', 'gslcblas'],
}

hdf5_dir = '/home/borg/smutch/cephfs/3rd_party/hdf5-1.8.17'
deps['hdf5'] = {
    'inclp' : os.path.join(hdf5_dir, 'include'),
    'libp'  : os.path.join(hdf5_dir, 'lib'),
    'lib'   : ['hdf5', 'hdf5_hl', 'z'],
}

gbp_dir = '/home/borg/smutch/cephfs/3rd_party/gbpCode'
deps['gbpCode'] = {
    'inclp' : os.path.join(gbp_dir, 'myInclude'),
    'libp'  : os.path.join(gbp_dir, 'myLib/mpi'),
    'lib'   : 'gbpLib',
}

fftw_dir = '/home/borg/smutch/cephfs/3rd_party/fftw-3.3.5'
deps['fftw'] = {
    'inclp' : os.path.join(fftw_dir, 'include'),
    'libp'  : os.path.join(fftw_dir, 'lib'),
    'lib'   : ['fftw3f_mpi', 'fftw3f'],
}

zlog_dir = '/home/borg/smutch/cephfs/3rd_party/zlog'
deps['zlog'] = {
    'inclp' : os.path.join(zlog_dir, 'include'),
    'libp'  : os.path.join(zlog_dir, 'lib'),
    'lib' : 'zlog',
}
