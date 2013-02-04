#!/usr/bin/env python

from fabricate import *
import sys
from glob import glob
import optparse

comiler = 'gcc'
sources = ['test_read', 'test_propogate', 'test_walk',]
flags = '-Wall -O2 -std=gnu99'.split()
# flags = '-Wall -O2'.split()

extra_opts = [optparse.make_option('-o','--options', type='string',
                                   help='Build options',
                                   default=''),]

if sys.platform == 'linux2':
    packages = {
        "cn_exceptions" : {
            "incpath" : "/home/smutch/code/chuck_norris_exceptions/include",
            "libpath" : "/home/smutch/code/chuck_norris_exceptions/lib",
        },
        "gbpLib" : {
            "incpath" : "/home/smutch/3rd_party/gbpCode/myInclude",
            "libpath" : "/home/smutch/3rd_party/gbpCode/myLib",
        },
    }
    ctags_exec = '/usr/bin/ctags'
else:
    packages = {
        "cn_exceptions" : {
            "incpath" : "/Users/smutch/Code/my_packages/chuck_norris_exceptions/include",
            "libpath" : "/Users/smutch/Code/my_packages/chuck_norris_exceptions/lib",
        },
        "gbpLib" : {
            "incpath" : "/Users/smutch/3rd_party/gbpCode/myInclude",
            "libpath" : "/Users/smutch/3rd_party/gbpCode/myLib",
        },
    }
    ctags_exec = '/usr/local/bin/ctags'

incpaths = ['-I'+packages[p]["incpath"] for p in packages.keys()]
libpaths = ['-L'+packages[p]["libpath"] for p in packages.keys()]
libs = ['-l'+p for p in packages.keys()]
    
def debug():
    flags.remove('-O2')
    flags.append('-DDEBUG -O0 -g'.split())
    sources = ['test_walk', 'test_read',]
    build(sources)

def test_walk():
    sources = ['test_walk',]
    build(sources)

def build(sources=sources):
    for source in sources:
        run(comiler, '-o', source, 'read_trees.c', source+'.c', flags,
            incpaths, libpaths, libs, main.options.options.split())

def clean():
    autoclean()

def ctags():
    run(ctags_exec, glob("*.[ch]"), "-R --c++-kinds=+p --fields=+iaS --extra=+q".split())

main(parallel_ok=True, extra_options=extra_opts, jobs=3)
