Installation
============

The TEDS package is written in Python and C++. All modules are available in Python and the L1A-L1B processor and instrument model are additionally also in C++:

.. csv-table::
   :align: center
   :header: Python, C++
   :widths: auto

   Geometry module
   Scene generator
   Instrument model, Instrument model
   L1A-L1B processor, L1A-L1B processor
   L1B-L2 processor
   L2-L4 processor

Installation of the Python modules means setting up a Python virtual environment using the pip package management system. The C++ code, on the other hand, needs to be compiled by the user after all the dependencies have been installed. This requires some knowledge of the build system as described in this chapter.

TEDS has only been tested on Linux operating systems. MacOS and Windows are not officially supported although in principle all prerequisites (libraries and compilers) also exist for those. In the following, we demonstrate how to run some of the commands on Ubuntu 24.04.


Python code
-----------

In order to ensure compatibility with the correct Python packages, it is important to work in a Python virtual environment. First install some prerequisites:

.. code-block:: bash

   sudo apt-get install python3.12-venv libpython3-dev libnetcdf-c++4-dev g++

Some of these are required for building the C++ extensions of TEDS. Now you can create a virtual environment with

.. code-block:: bash

   python -m venv venv  # Correct command might be python3

and activate it with

.. code-block:: bash

   source venv/bin/activate

You will notice that the shell's prompt has changed to remind you that you are in a virtual environment. Any packages installed with the Python :command:`pip` command are now part of the current project only. Correct versions of packages that are required for this project are listed in :file:`pyproject.toml`. Install them by issuing

.. code-block:: bash

   pip install --upgrade pip
   pip install --editable .

from the root source directory. The second command installs all dependencies found in :file:`pyproject.toml` and creates an *editable* build of TEDS suitable for development. Here the C++ environment is defined by the initial cache file that is in the root source directory. If you want to override that run

.. code-block:: bash

   INITIAL_CACHE_FILE=initial_cache.cmake pip install -editable .

If your working directory is the same as where the virtual environment is located, the TEDS package should automatically be found by Python scripts. If not, you might need to update your ``PYTHONPATH``.


C++ code
---------


Prerequisites
+++++++++++++

The instrument model and the L1A-L1B processor are dependent on the CMake build system generator. CMake is a set of tools for configuring, building, and testing software, and is released under the New BSD License. Building the processor is a two-step process: first run CMake to generate the build files, e.g. GNU Makefiles, and then issue the build command, e.g. make.

First step is to obtain CMake. On Ubuntu and similar distributions it can likely be installed by issuing

.. code-block:: bash

   sudo apt install cmake

Running CMake (the configuration step) creates a set of persistent variables which are contained in a file called :file:`CMakeCache.txt` in the build directory. These are referred to as cache variables and they are the user-configurable settings of the project. All the important decisions such as which compiler to use, which libraries to link against, etc., are stored as cache variables. There are several ways of setting the cache variables, one of which is to define them in a file that can be read by CMake. This is called the initial cache file, template of which are provided with the source code so you don't have to compose it from scratch.

The C++ code depends on the following libraries:

 * spdlog -- a popular C++ logging library
 * yaml-cpp -- a YAML parser
 * NetCDF -- self-describing data format library
 * pocketfft -- library for fast Fourier transforms
 * Eigen3 -- linear algebra library

You can install the first four with your Linux distribution's package manager. For example, on Ubuntu,

.. code-block:: bash

   sudo apt install libspdlog-dev libyaml-cpp-dev libnetcdf-c++4-dev libeigen3-dev

A copy of pocketfft is hosted at Bitbucket. You can clone it with

.. code-block:: bash

   git clone git@bitbucket.org:sron_earth/pocketfft.git

That said, you only need to ensure that NetCDF and Eigen are present. The rest, if not found, are downloaded and built automatically.

Both C++ codes depend on an OpenMP capable C++ compiler is required. Any recent version of the GNU C++ compiler :command:`g++` will do. If not already present, install by issuing

.. code-block:: bash

   sudo apt install g++


Configure and build
+++++++++++++++++++++

Most of the C++ code resides in the L1A-L1B processor and the instrument uses it as a dependency. A CMakeLists.txt found in the root source directory is a CMake script that binds them into a single project.

Start by navigating into the source directory and make a copy of the initial cache file:

.. code-block:: bash

   cd <teds>
   cp initial_cache.cmake initial_cache_local.cmake

where :file:`<teds>` denotes the root source directory of the TEDS project. Next, edit the initial cache file to reflect your environment, although the default values might already be fine (in which case there is no need to make a local copy of the file). When done editing, create a build directory and run CMake from that using the initial cache file:

.. code-block:: bash

   mkdir build && cd build
   cmake -C ../initial_cache_local.cmake ..

One can also build directly in the source directory but it is generally a good practice to do out-of-source builds and keep the source directory clean.

Note that editing the initial cache file has no effect after the first configuring! Instead, it is necessary to empty the build directory before running CMake again:

.. code-block:: bash

   rm -rf * # From the build directory
   cmake -C ../initial_cache_local.cmake ..

.. tip::

   Alternatively, if you want to keep the build directory intact while editing a CMake cache variable such as a compiler flag or a library to be linked against, you can use a graphical CMake front end or specify a given variable from the command line (the latter will not be demonstrated here). The two commonly used graphical front ends are the command line based :command:`ccmake` and the Qt-based :command:`cmake-gui`, obtained by issuing

   .. code-block:: bash

      sudo apt-get install cmake-curses-gui
      # or
      sudo apt-get install cmake-gui

   When using :command:`ccmake` issue

   .. code-block:: bash

      ccmake .

   from the build directory. Some CMake variables and options appear, most of which should be self-explanatory. A short help text to each variable is displayed at the bottom in a status bar. Pressing :kbd:`t` reveals all options. When done editing, press :kbd:`c` to reconfigure and :kbd:`g` to generate the native build files. Pay attention when :command:`ccmake` warns you that the cache variables have been reset. This will happen, e.g., when changing the compiler, and will necessitate the reconfiguring of some variables.

If CMake ran successfully the next step is to compile the executable. The default build system generated by CMake is GNU makefiles on Linux. Unless you are using a different build system, you can compile with

.. code-block:: bash

   make -j # or make -j VERBOSE=1 for more verbose output

If you are not sure which build system you are using, run

.. code-block:: bash

   cmake --build . # make is probably fine though

from the build directory. If successful, an executables called :file:`tango_l1b.x` and :file:`tango_im.x` are produced in the build directory.

The L1A-L1B processor can also be built independently because unlike the instrument model, it forms part of the operational processor. For that, navigate into its sources directory :file:`<teds>/teds/l1al1b` and follow the same steps as above. If all went well then only the :file:`tango_l1b.x` executable is produced.

.. tip::

   A different build system can be chosen by passing an argument to the CMake generator function. For instance, for using Ninja, use :command:`-G Ninja` during the initial configuring,

   .. code-block:: bash

      cmake -G Ninja -C <im>/initial_cache.cmake <im>

   The build command is then

   .. code-block:: bash

      ninja
      # or
      cmake --build .

CMake configuration variables
+++++++++++++++++++++++++++++++

See :file:`initial_cache.cmake` in the root directory of either C++ code for a list of configuration variables. You can copy and work with that file directly. There is thus no need to list them separately here.
