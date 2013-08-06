#!/usr/bin/env python

"""A collection of functions for doing common processing tasks on the galaxies
output by Meraxes."""

import numpy as np
from astropy import log

__author__ = "Simon Mutch"
__date__   = "2013/08/06"

def mass_function(mass, volume, **kwargs):
    """Generate a mass function.
    
    Args:
        mass   = an array of masses
        volume = volume of simulation cube / subset
        kwargs = passed to np.histogram call

    Returns:
        An array of [bin centers, mass function vals].
    """

    if "normed" in kwargs:
        kwargs["normed"] = False
        log.warn("Turned of normed kwarg in mass_function()")

    vals, edges = np.histogram(mass, **kwargs)
    width = edges[1]-edges[0]
    radius = width/2.0
    centers = edges[:-1]+radius
    vals = vals.astype(float) / (volume * width)

    return np.dstack((centers, vals)).squeeze()

