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

# Attempt to find the NetCDF library depending on which package
# manager is used
try:
    # Determine package name and location by searching for 'netcdf'
    res = subprocess.check_output(['dpkg', '-S', 'netcdf']).decode('utf-8')
    netcdf_header = list(
        filter(lambda x: x.endswith('/netcdf'), res.split()))[0]
    netcdf_library = list(
        filter(lambda x: x.endswith('/libnetcdf_c++4.so'), res.split()))[0]
    include_dirs.append(os.path.dirname(netcdf_header))
    library_dirs.append(os.path.dirname(netcdf_library))
except FileNotFoundError or subprocess.CalledProcessError:
    pass

extension_geo = Extension('teds.l1al1b.python.geolocation',
                          sources,
                          include_dirs=include_dirs,
                          library_dirs=library_dirs,
                          libraries=['netcdf_c++4'],
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
