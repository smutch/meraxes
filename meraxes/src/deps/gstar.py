# Gstar library dependencies

deps = {}

deps['exec'] = {
    'mpicc' : '/usr/local/x86_64/gnu/openmpi-1.6.1/bin/mpicc',
    'cc' : 'cc',
    'git' : '/usr/bin/git',
}

deps['gsl'] = {
    'inclp' :'/usr/local/x86_64/gnu/gsl-1.9/include',
    'libp'  :'/usr/local/x86_64/gnu/gsl-1.9/lib',
    'lib'   : ['gsl', 'gslcblas'],
}

deps['hdf5'] = {
    'inclp' :'/usr/local/x86_64/gnu/hdf5-1.8.14-openmpi/include',
    'libp'  :'/usr/local/x86_64/gnu/hdf5-1.8.14-openmpi/lib',
    'lib'   : ['hdf5', 'hdf5_hl', 'z'],
}

deps['gbpCode'] = {
    'inclp' :'/home/smutch/3rd_party/gbpCode/myInclude',
    'libp'  :'/home/smutch/3rd_party/gbpCode/myLib',
    'lib'   : 'gbpLib',
}

deps['fftw'] = {
    'inclp' :'/home/smutch/3rd_party/fftw-3.3.3/include',
    'libp'  :'/home/smutch/3rd_party/fftw-3.3.3/lib',
    'lib'   : ['fftw3f_mpi', 'fftw3f'],
}

deps['cmocka'] = {
    'inclp' : '/home/smutch/3rd_party/cmocka/include',
    'libp' : '/home/smutch/3rd_party/cmocka/lib',
    'lib' : 'cmocka',
}
