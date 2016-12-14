# COEPP library dependencies
import os

deps = {}

deps['exec'] = {
    'mpicc' : '/cvmfs/experimental.cloud.coepp.org.au/sw/sl6/gcc484/x86_64/openmpi/2.0.1_tcp/bin/mpicc',
    'cc' : '/usr/bin/gcc',
    'git' : '/usr/bin/git',
}

gsl_dir = '/cvmfs/experimental.cloud.coepp.org.au/sw/sl6/gcc484/x86_64/gsl/2.1'
deps['gsl'] = {
    'inclp' : os.path.join(gsl_dir, 'include'),
    'libp'  : os.path.join(gsl_dir, 'lib'),
    'lib'   : ['gsl', 'gslcblas'],
}

hdf5_dir = '/coepp/cephfs/dragons/smutch/3rd_party/hdf5-1.10.0-patch1'
deps['hdf5'] = {
    'inclp' : os.path.join(hdf5_dir, 'include'),
    'libp'  : os.path.join(hdf5_dir, 'lib'),
    'lib'   : ['hdf5', 'hdf5_hl', 'z'],
}

gbp_dir = '/coepp/cephfs/dragons/smutch/3rd_party/gbpCode'
deps['gbpCode'] = {
    'inclp' : os.path.join(gbp_dir, 'myInclude'),
    'libp'  : os.path.join(gbp_dir, 'myLib/mpi'),
    'lib'   : 'gbpLib',
}

fftw_dir = '/coepp/cephfs/dragons/smutch/3rd_party/fftw-3.3.5'
deps['fftw'] = {
    'inclp' : os.path.join(fftw_dir, 'include'),
    'libp'  : os.path.join(fftw_dir, 'lib'),
    'lib'   : ['fftw3f_mpi', 'fftw3f'],
}

zlog_dir = '/coepp/cephfs/dragons/smutch/3rd_party/zlog'
deps['zlog'] = {
    'inclp' : os.path.join(zlog_dir, 'include'),
    'libp'  : os.path.join(zlog_dir, 'lib'),
    'lib' : 'zlog',
}
