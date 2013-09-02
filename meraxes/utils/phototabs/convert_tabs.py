#!/usr/bin/env python

"""Convert ascii tables to hdf5 for use with Meraxes."""

import h5py as h5
import numpy as np
import pandas as pd
from glob import glob
import re

__author__ = "Simon Mutch"
__date__   = "2013/09/02"

# # conroy tables
# bands = ("U", "Bub", "Bbv", "V", "Rc", "Ic", "J", "H", "K", "u", "g", "r", "i", "z")
# for imf in ("salpeter", "kroupa", "chabrier"):
#     for mag_system in ("vega", "AB"):
#         base_dir = "conroy09/"+imf+"/"+mag_system+"/"
#         ages = pd.read_csv(base_dir+"conroy_agesMyr.dat", skiprows=1, names=('age',), header=None)

#         fnames = glob(base_dir+"mags*.c09")
#         metal_vals = [float(re.search('Z0\.[0-9]+', f).group()[1:]) for f in fnames]
#         phototabs = pd.DataFrame()
#         for val, f in zip(metal_vals, fnames):
#             temp = pd.read_csv(f, sep='    ', names=bands, header=None, )
#             temp["metallicity"] = pd.Series(np.repeat(val, len(temp.index)), name="metallicity")
#             temp = temp.join(ages)
#             phototabs = phototabs.append(temp)
#         metallicities = phototabs.groupby('metallicity')

#         with h5.File("../../input/photometric_tables/conroy09/"+imf+"/"+mag_system+".hdf5", "w") as fout:
#            for b in bands:
#                group_out = fout.create_group("%s"%b)
#                for name, g in metallicities:
#                    # print "%s::%s -> %.4f"%(imf, mag_system, name)
#                    group_out.create_dataset("%.04f"%name, data=g[b], dtype='f4')
#            fout.create_dataset("age", data=ages.values[:,0], dtype='f4')
#            fout["age"].attrs["units"] = "Myr" 

# original (Croton+ 2006) bc03 vega tables
bands = ["Bbv", "V", "Rc", "I", "K"]
base_dir = "bc03/salpeter/vega/"
fnames = glob(base_dir+"*.bc03")
metal_vals = [float(re.search('Z0\.[0-9]+', f).group()[1:]) for f in fnames]
phototabs = pd.DataFrame()
for val, f in zip(metal_vals, fnames):
    temp = pd.read_csv(f, sep=' +', header=None)
    ages = temp[6].copy()
    temp = temp.drop(0, axis=1)
    temp.columns = bands+["age",]
    temp["metallicity"] = pd.Series(np.repeat(val, len(temp.index)), name="metallicity")
    phototabs = phototabs.append(temp)
metallicities = phototabs.groupby('metallicity')


with h5.File("../../input/photometric_tables/bc03/salpeter/vega.hdf5", "w") as fout:
   for b in bands:
       group_out = fout.create_group("%s"%b)
       for name, g in metallicities:
           group_out.create_dataset("%.04f"%name, data=g[b], dtype='f4')
   fout.create_dataset("age", data=ages.values, dtype='f4')
   fout["age"].attrs["units"] = "Myr" 

