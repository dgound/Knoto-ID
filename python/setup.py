import sys
import os
import subprocess
from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pathlib import Path

def find_boost_dirs():
    """Try to locate Boost include and lib directories automatically.

    Detection order: explicit environment variables (BOOST_ROOT, or
    BOOST_INCLUDEDIR / BOOST_LIBRARYDIR), then a conda environment, then
    Homebrew (macOS), then pkg-config, then common system paths. Set
    BOOST_ROOT (or BOOST_INCLUDEDIR / BOOST_LIBRARYDIR) to point at a Boost
    installation in a non-standard location, e.g. on a server or an HPC
    module.
    """
    boost_include = None
    boost_lib = None

    # 1. Explicit environment variables (highest priority).
    include_env = os.getenv('BOOST_INCLUDEDIR')
    lib_env = os.getenv('BOOST_LIBRARYDIR')
    boost_root = os.getenv('BOOST_ROOT') or os.getenv('BOOSTROOT')
    if include_env:
        boost_include = include_env
    elif boost_root:
        if (Path(boost_root) / "include" / "boost").exists():
            boost_include = str(Path(boost_root) / "include")
        elif (Path(boost_root) / "boost").exists():
            boost_include = boost_root
    if lib_env:
        boost_lib = lib_env
    elif boost_root:
        for cand in ("lib", "lib64", "stage/lib"):
            if (Path(boost_root) / cand).exists():
                boost_lib = str(Path(boost_root) / cand)
                break

    # 2. conda environment.
    conda_prefix = os.getenv('CONDA_PREFIX')
    if conda_prefix:
        conda_include = Path(conda_prefix) / "include"
        conda_lib = Path(conda_prefix) / "lib"
        if not boost_include and (conda_include / "boost").exists():
            boost_include = str(conda_include)
        if not boost_lib and conda_lib.exists():
            boost_lib = str(conda_lib)

    # 3. macOS: Boost installed via Homebrew.
    if sys.platform == "darwin" and not boost_include:
        try:
            homebrew_prefix = subprocess.check_output(
                ['brew', '--prefix'], universal_newlines=True
            ).strip()
            hb_include = Path(homebrew_prefix) / "include"
            hb_lib = Path(homebrew_prefix) / "lib"
            if (hb_include / "boost").exists():
                boost_include = str(hb_include)
                if not boost_lib and hb_lib.exists():
                    boost_lib = str(hb_lib)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

    # 4. pkg-config (only fills in what is still missing).
    if not boost_include or not boost_lib:
        try:
            inc = subprocess.check_output(
                ['pkg-config', '--cflags-only-I', 'boost'],
                universal_newlines=True
            ).strip().replace("-I", "")
            lib = subprocess.check_output(
                ['pkg-config', '--libs-only-L', 'boost'],
                universal_newlines=True
            ).strip().replace("-L", "")
            if inc and not boost_include:
                boost_include = inc
            if lib and not boost_lib:
                boost_lib = lib
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

    # 5. Common system paths. Boost headers live in /usr/include/boost, so the
    #    include directory is /usr/include (the parent), not /usr/include/boost.
    if not boost_include:
        boost_include = "/usr/include"
    if not boost_lib:
        boost_lib = "/usr/lib"

    return boost_include, boost_lib


# Automatically find Boost directories
boost_include, boost_lib = find_boost_dirs()

# Include directory locations
include_dirs = [
    ".",
    "KnotoID/cpp_bindings/",
    "../src",
    boost_include
]

# Library directories
library_dirs = [boost_lib]

# Define the Pybind11 extension
ext_modules = [
    Pybind11Extension(
        "knotoID_cpp",
        [
            "KnotoID/cpp_bindings/bindings.cpp",
            "KnotoID/cpp_bindings/KnottedCore.cpp",
            "../src/PlanarDiagram.cpp",
            "../src/PlanarGraph.cpp",
            "../src/Polygon.cpp",
            "../src/Polynomial.cpp",
            "../src/PolynomialInvariant.cpp",
            "../src/DoubleBranchedCover.cpp",
            "../src/Random.cpp"
        ],
        include_dirs=include_dirs,
        library_dirs=library_dirs,
    ),
]

# Generic setup for cross-platform compatibility
setup(
    name="KnotoID",
    version="1.4.0",
    author="Julien Dorier, Dimos Goundaroulis",
    description="Python bindings for Knoto-ID",
    packages=find_packages(),
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    install_requires=[
        "numpy",
        "pandas",
        "matplotlib",
        "plotnine",  # ggplot2-compatible rendering for plot_knotted_core
    ],
    package_data={
        'KnotoID': ['data/*.txt'],  # Include all txt files from data directory
    }
)
