Meraxes
=======

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit)

The Meraxes semi-analytic galaxy formation model.

> [!WARNING]
> This code is in pre-release phase and contains unpublished physical recipes and code implementations.  Please do not distribute this code without permission or copy any portions.  All rights reserved!

> [!TIP]
> See `BUILD.md` for detailed build instructions.


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

If you use any of the following features, please also cite the relevant papers:

AGN model
: [Qin et al. 2017 (arXiv:1703.04895)](https://ui.adsabs.harvard.edu/#abs/2017arXiv170304895Q/abstract)

Minihalo model
: [Ventura et al. Subm. (arXiv:2401.07396)](https://ui.adsabs.harvard.edu/abs/2024arXiv240107396V/abstract)

Contributors
------------

* **Simon Mutch** (Lead author; *The University of Melbourne*)
* Paul Geil (*The University of Melbourne*)
* Yuxiang Qin (*The University of Melbourne*)
* Hansik Kim (*The University of Melbourne*)
* Greg Poole (*Astronomy Data and Compute Services*)
* Yisheng Qiu (*The University of Melbourne*)
* Brad Greig (*The University of Melbourne*)
* Emanuele Maria Ventura (*The University of Melbourne*)


How to contribute
=================

1. Fork the [official repository](https://github.com/smutch/meraxes).
2. Make your changes in a new, appropriately named branch (e.g. `mutch_starformation` or `fix_progenitor_order`).
    * All git commit messages should loosely follow the [standard format outlined here](http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html).
    * Ensure your new code is commented enough so that others can understand what it's doing.
3. Issue a [pull request](https://help.github.com/en/articles/about-pull-requests) to the official repository.
