# Instructions for Setuptools for building the C++ component TEDS

from setuptools import Extension
from setuptools import setup
from setuptools.command.build_ext import build_ext
import os
import platform
import site
import subprocess

if platform.system() == 'Windows':
    setup(ext_modules=[])
    exit(0)


class CMakeExtension(Extension):
    def __init__(self, name):
        super().__init__(name,
                         sources=[],
                         py_limited_api=True)
        self.source_dir = os.path.abspath('')


class CMakeBuild(build_ext):
    def build_cmake(self, ext):
        if 'INITIAL_CACHE_FILE' in os.environ:
            initial_cache_file = os.environ['INITIAL_CACHE_FILE']
        else:
            initial_cache_file = ext.source_dir + '/initial_cache.cmake'
        lib_dest_dir = site.getsitepackages()[0] + '/teds_cpp'
        cmake_args = [f'-C {initial_cache_file}',
                      f'-D CMAKE_LIBRARY_OUTPUT_DIRECTORY={lib_dest_dir}']
        if subprocess.call(['which', 'ninja']) == 0:
            cmake_args += ['-G Ninja']
        build_args = ["--config", "Release"]
        self.build_lib = self.build_lib
        os.makedirs(self.build_lib, exist_ok=True)
        subprocess.check_call(
            ["cmake", ext.source_dir] + cmake_args, cwd=self.build_lib)
        subprocess.check_call(
            ["cmake", "--build", "."] + build_args, cwd=self.build_lib)

    def run(self):
        for ext in self.extensions:
            self.build_cmake(ext)


setup(ext_modules=[CMakeExtension('cpp.bindings')],
      cmdclass={'build_ext': CMakeBuild})
