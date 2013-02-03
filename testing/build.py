from fabricate import *

comiler = 'mpicc'
sources = ['test_read', 'test_propogate', 'test_walk',]
flags = '-Wall -O2 -std=gnu99 -DDEBUG'.split()
# flags = '-Wall -O2'.split()
ctags_exec = '/usr/local/bin/ctags'

packages = {
    "cn_exceptions" : {
        "incpath" : "/Users/smutch/Code/my_packages/chuck_norris_exceptions/include",
        "libpath" : "/Users/smutch/Code/my_packages/chuck_norris_exceptions/lib",
    },
    "gbpLib" : {
        "incpath" : "/Users/smutch/3rd_party/gbpCode/myInclude",
        "libpath" : "/Users/smutch/3rd_party/gbpCode/myLib/mpi",
    },
}

incpaths = ['-I'+packages[p]["incpath"] for p in packages.keys()]
libpaths = ['-L'+packages[p]["libpath"] for p in packages.keys()]
libs = ['-l'+p for p in packages.keys()]

def debug():
    cflags.append('-g')
    cflags.remove('-O2')
    cflags.append('-O0')
    build()

def test_walk():
    sources = ['test_walk',]
    build(sources)

def build(sources=sources):
    for source in sources:
        run(comiler, '-o', source, 
            'read_trees.c', source+'.c',
            flags, incpaths, libpaths, libs)

def clean():
    autoclean()

def ctags():
    run(ctags_exec+" -R --c++-kinds=+p --fields=+iaS --extra=+q")

main(parallel_ok=True)
