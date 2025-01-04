# Instructions for Setuptools for building the Python geolocation
# interface with a Numpy dependency.

from setuptools import Extension
from setuptools import setup
import numpy

sources = ['teds/l1al1b/tango_l1b/algorithm.cpp',
           'teds/l1al1b/tango_l1b/quaternion.cpp',
           'teds/l1al1b/tango_l1b/solar_model.cpp',
           'teds/l1al1b/tango_l1b/dem.cpp',
           'teds/l1al1b/tango_l1b/geometry.cpp',
           'teds/l1al1b/python/geolocation.cpp']
extension = Extension('teds.l1al1b.python.geolocation',
                      sources,
                      include_dirs=[numpy.get_include()],
                      libraries=['hdf5'],
                      extra_compile_args=['-std=c++20'])
setup(ext_modules=[extension])
