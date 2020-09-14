Meraxes
=======

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

The Meraxes semi-analytic galaxy formation model.

**Author**: [Simon Mutch](http://www.ph.unimelb.edu.au/~smutch/index.html) (*The University of Melbourne*)

**Note: This code is in pre-release phase and contains unpublished physical recipes and code implementations.  Please do not distribute this code without permission or copy any portions.  All rights reserved!**

Requirements
------------

* [GSL](https://www.gnu.org/software/gsl/)
* [FFTW3](http://www.fftw.org) built with `--enable-float --enable-mpi` ([see installation docs](http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix)).
* [MLOG](https://github.com/smutch/mlog) (This is packaged with Meraxes as a git submodule and will be built automatically.)
* [HDF5](https://support.hdfgroup.org/HDF5/) built with MPI bindings.


Optional
--------

* A NVIDIA GPU and [CUDA](https://developer.nvidia.com/cuda-zone) to accelerate the reionisation calculation


Attribution
-----------

If you use this code or its results in any publications, please cite [Mutch et al. 2016b](http://adsabs.harvard.edu/abs/2016MNRAS.462..250M) (MNRAS, 462, 250):

```
@ARTICLE{2016MNRAS.462..250M,
   author = {{Mutch}, S.~J. and {Geil}, P.~M. and {Poole}, G.~B. and {Angel}, P.~W. and 
	{Duffy}, A.~R. and {Mesinger}, A. and {Wyithe}, J.~S.~B.},
    title = "{Dark-ages reionization and galaxy formation simulation - III. Modelling galaxy formation and the epoch of reionization}",
  journal = {\mnras},
archivePrefix = "arXiv",
   eprint = {1512.00562},
 keywords = {galaxies: formation, galaxies: high redshift, dark ages, reionization, first stars},
     year = 2016,
    month = oct,
   volume = 462,
    pages = {250-276},
      doi = {10.1093/mnras/stw1506},
   adsurl = {http://adsabs.harvard.edu/abs/2016MNRAS.462..250M},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

For the AGN model, please cite [Qin et al. 2017 (arXiv:1703.04895)](https://ui.adsabs.harvard.edu/#abs/2017arXiv170304895Q/abstract).


Contributors
------------

* Paul Geil (*The University of Melbourne*)
* Yuxiang Qin (*The University of Melbourne*)
* Hansik Kim (*The University of Melbourne*)
* Greg Poole (*Astronomy Data and Compute Services*)
* Yisheng Qiu (*The University of Melbourne*)
* Brad Greig (*The University of Melbourne*)


How to contribute
=================

1. Fork the [official repository](https://bitbucket.org/dragons-astro/meraxes).
2. Make your changes in a new, appropriately named branch (e.g. `mutch_starformation` or `fix_progenitor_order`).
    * All git commit messages should loosely follow the [standard format outlined here](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html).
    * Ensure your new code is commented enough so that others can understand what it's doing.
3. Issue a [pull request](https://help.github.com/en/articles/about-pull-requests) to the official repository.
