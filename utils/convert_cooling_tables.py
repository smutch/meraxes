#!/usr/bin/env python

"""Convert the ascii cooling function tables of Sutherland & Dopita (1993) to
hdf5."""

import pandas as pd
import h5py as h5
from os import path
import re

__author__ = "Simon Mutch"
__date__ = "2014/03/07"

regex = re.compile('CoolFunctions/(?P<m>.*).cie')


def read_ascii_cooling_file(fin):
    """Read in the input ascii file"""

    print "reading file", fin

    column_names = ("log(T)", "ne", "nH", "nt", "log(lambda_net)",
                    "log(lambda_norm)", "log(U)", "log(taucool)", "P12",
                    "rho24", "Ci", "mubar")

    data = pd.read_table(fin, delim_whitespace=True, skiprows=4,
                         header=None, names=column_names, dtype='f4')

    data = data.dropna()

    # If we don't have the first two rows (T=1e4) then just copy the first
    # present one
    if data.shape[0] < 90:
        first_val = data.ix[0].copy()
        first_val["log(T)"] = 4.05
        data = pd.concat((pd.DataFrame(first_val).T, data), axis=0)
        first_val["log(T)"] = 4.0
        data = pd.concat((pd.DataFrame(first_val).T, data), axis=0)

    return data


def save_data_to_hdf5(fout, data, fin):
    """Save the data to the output hdf5 file"""

    metallicity = regex.search(fin).group('m')
    print "creating group", metallicity
    group = fout.create_group(metallicity)

    for col in data.columns:
        group.create_dataset(col, data=data[col])


def append_metadata(fout):

    meta = {
        "mzero": "cie cooling for primordial hydrogen/helium mix"
        "log(T) = 4-8.5",
        "m-00": "cie cooling for solar abundances mix log(T) = 4-8.5",
        "m-05": "[Fe/H] = -0.5, solar/primordial average ratios",
        "m-10": "[Fe/H] = -1.0, primordial ratios (ie enhanced oxygen)",
        "m-15": "[Fe/H] = -1.5, primordial ratios",
        "m-20": "[Fe/H] = -2.0, primordial ratios",
        "m-30": "[Fe/H] = -3.0, primordial ratios",
        "m+05": "[Fe/H] = +0.5, solar ratios log(T) = 4.1-8.5 (due to charge"
        "exchange problems at log(T) = 4.0)",
    }

    columns_description = """log(T): log temperature in K

ne, nH, nt : number densities, electrons, hydrogen and total
ion in cm^-3.

log(lambda net) log(lambda norm) : log on the net cooling function
and the normalised cooling function.  lambda norm = lambda net / (ne nt).
lambda net in ergs cm^-3 s^-1,  lambda net in ergs cm^3 s^-1.  While the
densities are kept less than about p/k 10^8 both isobaric and isochoric
curves can be constructed from the normalised function using an appropriate
density function.  The non-equilibrium net cooling function is from the
isobaric model used to calculate the curves.  In the CIE curves the net
function is for the isochoric cie model.

log(U):  U = 3/2 N kT, N = ne + nt the total internal energy. ergs cm^-3

log(taucool):  The normalized cooling timescale Nr*( U/(lambda net))
# Nr = (ne nt)/(ne+nt).  s cm^-3

P12:  Gas pressure NkT. ergs ergs cm^-3 times 10^12

rho24:  Density g  cm^-3 times 10^24

Ci: The isothermal sound speed kms s^-1

mubar: mean molecular weight grams times 10^24"""

    for k, v in meta.iteritems():
        fout[k].attrs["description"] = v

    fout.attrs["columns"] = columns_description


if __name__ == '__main__':

    input_dir = "/home/smutch/models/CoolFunctions"
    input_flist = ("m+05.cie", "m-00.cie", "m-05.cie", "m-10.cie", "m-15.cie",
                   "m-20.cie", "m-30.cie", "mzero.cie")

    # open the output file
    fout = h5.File("../input/cooling_functions/SD93.hdf5")

    # write a brief description
    description = """Collisional Ionisation Equilibrium cooling functions taken
    from Sutherland & Dopita (1993) and converted from the ascii input used by
    the Croton+ 2006 semi-analytic model."""
    fout.attrs["description"] = description

    # loop through each input file and save it in the hdf5 file
    for fin in [path.join(input_dir, f) for f in input_flist]:
        data = read_ascii_cooling_file(fin)
        save_data_to_hdf5(fout, data, fin)

    # append metadata
    append_metadata(fout)

    # close the output file
    fout.close()
