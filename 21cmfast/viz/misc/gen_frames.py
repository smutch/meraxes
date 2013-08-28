#!/usr/bin/env python

from glob import glob
from littleworkers import Pool

__author__ = "Simon Mutch"
__date__   = "2013/06/07"

xH_DIR = "/lustre/projects/p004_swin/smutch/21cmfast_runs/Tiamat/resampled/Boxes" 
SNAPLIST = range(40,100)
SLICE = ":,:,128"
SLICE_DEPTH = "1"
SIM = "Tiamat_full"

xHfiles = sorted(glob(xH_DIR+"/xH*"))[::-1]

commands = ["python show_slice_with_density.py "+fname+" frames/frame_{:03d}.png "
            .format(ii)+SLICE+" "+SLICE_DEPTH+" "+SIM+" "+str(SNAPLIST[ii])+" --color_bar"
            for ii, fname in enumerate(xHfiles)]

lil = Pool(workers=3)

# bar = ProgressBar(len(commands))
# def update(proc):
#     bar.update()
# lil.run(commands, callback=update)
# bar.__exit__(None,None,None)

lil.run(commands, progressbar=True)
