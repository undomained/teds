# Instructions for Setuptools for building TEDS C++ extensions

from setuptools import Extension
from setuptools import setup
import numpy
import os
import subprocess

sources = ['teds/l1al1b/tango_l1b/algorithm.cpp',
           'teds/l1al1b/tango_l1b/quaternion.cpp',
           'teds/l1al1b/tango_l1b/solar_model.cpp',
           'teds/l1al1b/tango_l1b/dem.cpp',
           'teds/l1al1b/tango_l1b/geometry.cpp',
           'teds/l1al1b/python/geolocation.cpp']
include_dirs = [numpy.get_include()]
library_dirs = []
compiler_flags = ['-O3', '-std=c++20', '-fopenmp']

# Attempt to find the HDF5 library depending on which package manager is used
try:
    # Determine package name and location by searching for hdf5.h
    res = subprocess.check_output(['dpkg', '-S', 'hdf5.h']).decode('utf-8')
    package_name, hdf5_header = res.split()
    package_name = package_name[:-1]
    include_dirs.append(os.path.dirname(hdf5_header))
    # Find which directory contains the HDF5 library
    res = subprocess.check_output(['dpkg', '-L', package_name]).decode('utf-8')
    for line in res.split():
        if 'libhdf5.so' in line:
            library_dirs.append(os.path.dirname(line))
            break
except FileNotFoundError or subprocess.CalledProcessError:
    pass

extension_geo = Extension('teds.l1al1b.python.geolocation',
                          sources,
                          include_dirs=include_dirs,
                          library_dirs=library_dirs,
                          libraries=['hdf5'],
                          extra_compile_args=compiler_flags,
                          extra_link_args=['-lgomp'],
                          py_limited_api=True)

sources = ['teds/lib/algorithms.cpp']
extension_alg = Extension('teds.lib.algorithms',
                          sources,
                          include_dirs=include_dirs,
                          extra_compile_args=compiler_flags,
                          extra_link_args=['-lgomp'],
                          py_limited_api=True)

setup(ext_modules=[extension_geo, extension_alg])
