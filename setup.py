#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os.path
import pybind11 as py11
import configparser
from setuptools import find_packages
from packaging.version import LegacyVersion
from skbuild.exceptions import SKBuildError
from skbuild.cmaker import get_cmake_version

try:
    from skbuild import setup
except ImportError:
    sys.stderr.write(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in"
        " pyproject.toml yourself"
    )
    raise

# This is needed for versioneer to be importable when building with PEP 517.
# See <https://github.com/warner/python-versioneer/issues/193> and links
# therein for more information.
source_folder = os.path.dirname(__file__)
sys.path.append(source_folder)

parser = configparser.ConfigParser()
parser.read("setup.cfg")
cmake_args = ["-Dpybind11_DIR:PATH={}".format(py11.get_cmake_dir())]

_version = None
if ("cmake_config" in parser) and ("akantu_dir" in parser["cmake_config"]):
    sys.path.append(parser["cmake_config"]["akantu_dir"])
    try:
        import akantu_version

        _version = akantu_version.get_version()
    except ImportError:
        pass

try:
    import versioneer

    if not _version:
        _version = versioneer.get_version()
    setup_kw = {
        "version": _version,
        "cmdclass": versioneer.get_cmdclass(),
    }
    cmake_args.append("-DAKANTU_VERSION={}".format(_version))
except ImportError:
    # see https://github.com/warner/python-versioneer/issues/192
    print("WARNING: failed to import versioneer," " falling back to no version for now")
    setup_kw = {}


if "cmake_config" in parser:
    for k, v in parser["cmake_config"].items():
        k = k.upper()
        cmake_args.append("-D{}:BOOL={}".format(k, v))

akantu_libs = []

if "CI_AKANTU_INSTALL_PREFIX" in os.environ:
    akantu_libs.extend(
        [
            # paths comming from the manylinux install via gitlab-ci
            "/softs/view/lib/*",
            "/softs/view/lib64/*",
            os.path.join(os.environ["CI_AKANTU_INSTALL_PREFIX"], "lib64/*"),
            os.path.join(os.environ["CI_AKANTU_INSTALL_PREFIX"], "lib/*"),
        ]
    )
    cmake_args.extend(
        [
            "-DAKANTU_BYPASS_AKANTU_TARGET:BOOL=ON",
            "-DAkantu_DIR:PATH={}".format(
                os.path.join(
                    os.environ["CI_AKANTU_INSTALL_PREFIX"], "lib", "cmake", "Akantu"
                )
            ),
        ]
    )

# Add CMake as a build requirement if cmake is not installed or is too low a
# version
setup_requires = []
try:
    if LegacyVersion(get_cmake_version()) < LegacyVersion("3.4"):
        setup_requires.append("cmake")
except SKBuildError:
    setup_requires.append("cmake")


with open(os.path.join(source_folder, "README.md"), "r") as fh:
    long_description = fh.read()

setup(
    name="akantu",
    url="https://akantu.ch",
    author="Nicolas Richart",
    author_email="nicolas.richart@epfl.ch",
    description="Akantu: Swiss-Made Open-Source Finite-Element Library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    platforms="",
    license="L-GPLv3",
    license_files=["COPYING", "COPYING.lesser"],
    project_urls={
        "Bug Tracker": "https://github.com/akantu/akantu/issues",
    },
    setup_requires=setup_requires,
    install_requires=["numpy", "scipy"],
    package_data={"AkantuLibs": akantu_libs},
    packages=find_packages(where="python"),
    package_dir={"": "python"},
    include_package_data=False,
    cmake_args=cmake_args,
    cmake_languages=["CXX"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Natural Language :: English",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: C++",
        "Programming Language :: Python",
        "Topic :: Education",
        "Topic :: Scientific/Engineering",
    ],
    **setup_kw
)
