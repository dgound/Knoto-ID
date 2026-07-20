import sys
import os
import subprocess
from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pathlib import Path

def find_boost_dirs():
    """Try to locate Boost include and lib directories automatically."""
    boost_include = None
    boost_lib = None

    # Check if we are inside a conda environment
    conda_prefix = os.getenv('CONDA_PREFIX')
    if conda_prefix:
        conda_include = Path(conda_prefix) / "include"
        conda_lib = Path(conda_prefix) / "lib"
        
        if (conda_include / "boost").exists():
            boost_include = str(conda_include)
        if conda_lib.exists():
            boost_lib = str(conda_lib)

    # macOS specific handling for Boost installed via Homebrew
    if sys.platform == "darwin" and not boost_include:
        try:
            # Default Homebrew Boost installation path
            homebrew_prefix = subprocess.check_output(
                ['brew', '--prefix'], universal_newlines=True
            ).strip()
            boost_include = Path(homebrew_prefix) / "include"
            boost_lib = Path(homebrew_prefix) / "lib"

            if not (boost_include / "boost").exists():
                boost_include = None  # Boost is not installed
            if not boost_lib.exists():
                boost_lib = None
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

    # Fallback: Try to find Boost using pkg-config (if installed)
    if not boost_include or not boost_lib:
        try:
            boost_include = subprocess.check_output(
                ['pkg-config', '--cflags-only-I', 'boost'],
                universal_newlines=True
            ).strip().replace("-I", "")
            boost_lib = subprocess.check_output(
                ['pkg-config', '--libs-only-L', 'boost'],
                universal_newlines=True
            ).strip().replace("-L", "")
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass

    # Further fallback: Default to common system paths (Linux)
    if not boost_include:
        boost_include = "/usr/include/boost"
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
