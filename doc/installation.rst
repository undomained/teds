Installation
============

The TEDS package is written in Python and C++. The instrument model and L1A-L1B processor are in C++ and everything else in Python:

.. csv-table::
   :align: center
   :header: Python, C++
   :widths: auto

   Geometry module, Instrument model
   Scene generator, L1A-L1B processor
   L1B-L2 processor
   L2-L4 processor

Installation of the Python modules mainly refers to setting up a Python virtual environment using the pip package management system. The C++ code, once all the dependencies have been installed, needs to be compiled by the user which requires knowledge of the build system, described in this chapter.

We have tested only TEDS on Linux operating systems. MacOS and Windows are not officially supported but in principle all prerequisites (relevant libraries and compilers) also exist for those. In the following, we demonstrate how to run some of the commands only on Ubuntu 24.04.


Python code
-----------

Set up a virtual environment:

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate
   pip install --upgrade pip
   pip install -r teds/requirements


C++ code
---------


Prerequisites
+++++++++++++

The instrument model and the L1A-L1B processor are dependent on the CMake build system generator. CMake is a set of tools for configuring, building, and testing software, and is released under the New BSD License. Building the processor is a two-step process: first run CMake to generate the build files, e.g. GNU Makefiles, and then issue the build command, e.g. make.

First step is to obtain CMake. On Ubuntu and similar distributions it can likely be installed by issuing

.. code-block:: bash

   sudo apt install cmake

Running CMake (the configuration step) creates a set of persistent variables which are contained in a file called ``CMakeCache.txt`` in the build directory. These are referred to as cache variables and they are the user-configurable settings of the project. All the important decisions such as which compiler to use, which libraries to link against, etc., are stored as cache variables. There are several ways of setting the cache variables, one of which is to define them in a file that can be read by CMake. This is called the initial cache file, template of which are provided with the source code so you don't have to compose it from scratch.

Both C++ codes depend on the following libraries:

 * spdlog -- a popular C++ logging library
 * yaml-cpp -- a YAML parser
 * pocketfft -- library for fast Fourier transforms

You can install the first two with your Linux distribution's package manager. For example, on Ubuntu,

.. code-block:: bash

   apt-get install libspdlog-dev libyaml-cpp-dev

A copy of ``pocketfft`` is hosted at Bitbucket. You can clone it with

.. code-block:: bash

   git clone git@bitbucket.org:sron_earth/pocketfft.git

Also, an OpenMP capable C++ compiler is required. Any recent version of the GNU C++ compiler ``g++`` will do.


Configure and build
+++++++++++++++++++++

Most of the C++ code resides in the L1A-L1B processor and the instrument uses it as a dependency. Thus, the L1A-L1B processor needs to be built first.

Start by navigating into the source directory:

.. code-block:: bash

   cd <teds>/teds/L1AL1B/tango_l1b
   cp initial_cache.cmake.example initial_cache.cmake

where ``<teds>`` denotes the root source directory of the TEDS project. Next, edit the initial cache file to reflect your environment, although the default values might already be fine. When done editing, create a build directory and run CMake from that using the initial cache file:

.. code-block:: bash

   mkdir build && cd build
   cmake -C <teds>/teds/L1AL1B/tango_l1b/initial_cache.cmake <teds>/teds/L1AL1B/tango_l1b

One can also build directly in the source directory but it is generally a good habit to do out-of-source builds and keep the source directory clean.

Note that editing the initial cache file has no effect after the first configuring! Instead, it is necessary to empty the build directory before running CMake again:

.. code-block:: bash

   rm -rf * # From the build directory
   cmake -C <teds>/teds/L1AL1B/tango_l1b/initial_cache.cmake

.. tip::

   Alternatively, if you want to keep the build directory intact while editing a CMake cache variable such as a compiler flag or a library to be linked against, you can use a graphical CMake front end or specify that variable from the command line (the latter will not be demonstrated here). The two commonly used graphical front ends are the command line based `ccmake` and the Qt-based `cmake-gui`, obtained by issuing

   .. code-block:: bash

      sudo apt-get install cmake-curses-gui
      # or
      sudo apt-get install cmake-gui

   When using `ccmake` issue

   .. code-block:: bash

      ccmake .

   from the build directory. Some CMake variables and options appear, most of which should be self-explanatory. A short help text to each variable is displayed at the bottom in a status bar. Pressing `t` reveals all options. When done editing, press `c` to reconfigure and `g` to generate the native build files. Pay attention when `ccmake` warns you that the cache variables have been reset. This will happen, e.g., when changing the compiler, and will necessitate the reconfiguring of some variables.

If CMake ran successfully the next step is to compile the executable. The default build system generated by CMake is GNU makefiles on Linux. Unless you are using a different build system, you can compile with

.. code-block:: bash

   make -j # or make -j VERBOSE=1 for more verbose output

If you are not sure which build system you are using, run

.. code-block:: bash

   cmake --build . # make is probably fine though

from the build directory. If successful, an executable called `spexone_cal` is produced in the build directory.

.. tip::

   A different build system can be chosen by passing an argument to the CMake generator function. For instance, for using Ninja, use `-G Ninja` during the initial configuring,

   .. code-block:: bash

      cmake -G Ninja -C <spexone_cal>/initial_cache.cmake <spexone_cal>

   The build command is then

   .. code-block:: bash

      ninja
      # or
      cmake --build .

CMake configuration variables
-------------------------------

See `initial_cache.cmake.example` in the root directory for a list of configuration variables. You can copy and work with that file directly. There is thus no need to list them separately here.
