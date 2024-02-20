# Copyright 2013-2021 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack.build_systems.meson import MesonBuilder

from spack import *

MesonBuilder._std_args_orig = MesonBuilder.std_args


def meson_args(pkg):
    args = MesonBuilder._std_args_orig(pkg)
    args.remove("-Dwrap_mode=nodownload")
    args.append("-Dwrap_mode=default")
    return args

MesonBuilder.std_args = meson_args


class CriterionGit(MesonPackage):
    """A dead-simple, yet extensible, C and C++ unit testing framework."""

    homepage = "https://github.com/Snaipe/Criterion/tree/master"
    git = "git@github.com:Snaipe/Criterion.git"

    maintainers = [
        "smutch",
    ]

    version("2.4.2", tag="v2.4.2", submodules=True)

    # depends_on("gettext")
    depends_on("meson@0.55.0:", type="build")
    depends_on("ninja", type="build")
    depends_on("cmake", type="build")
    depends_on("pkg-config", type="build")
    depends_on("libffi", type="build")
    depends_on("libgit2", type="build")
