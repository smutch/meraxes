# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class CriterionGit(CMakePackage):
    """A dead-simple, yet extensible, C and C++ unit testing framework."""

    homepage = "https://github.com/Snaipe/Criterion/tree/master"
    git      = "git@github.com:Snaipe/Criterion.git"

    maintainers = ['smutch',]

    version('2.3.3', tag='v2.3.3', submodules=True)
