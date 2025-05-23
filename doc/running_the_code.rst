Running the code
====================


General remarks
---------------

TEDS has a modular design where each module focuses on one specific task in the end-to-end chain and has one of more entry points. This chapter describes how to run each module individually and provides an example script of a full-chain simulation. Ultimately, it's up to the user to write a driver script which contains instructions for how to execute all or some parts of the chain.


YAML configuration files
++++++++++++++++++++++++

Each TEDS module takes one argument which is a reference to a configuration file written in YAML. We chose YAML because it is versatile, human-readable, and widely used and supported. It uses Python-style indentation to indicate nesting and as basic types supports scalars, lists, and maps (also called dictionaries or hashes). Furthermore, it is a strict superset of JSON and thus compatible with it.

A quick way to understand the YAML format is run the L1A-L1B C++ code without any arguments:

.. code-block:: bash

   <build>/tango_l1b.x

The output contains keywords showing their type, default value, and description. A YAML file may contain several sections. For instance, a ``detector`` section contains a set of parameters specifically targeting the detector.

The basic element of YAML is a key-value pair where the value is a scalar, list, or a map. Here is example YAML where the element type is noted as a comment::

  io: # starts a map
    l1a_files: # key followed by a list
      - l1a_1.nc # scalar
      - l1a_2.nc # scalar
  cal_level: l1a # key-value

There are many checks in the code which notify you of errors and inconsistencies in the configuration file. If you use a non-existent keyword the processor will issue a warning but will continue with the calculation.

Full specification of YAML and how to use it is found here: https://yaml.org/spec/1.2.2/.


NetCDF
++++++++

It is an objective of TEDS to only work with NetCDF files, either as input or output of any module. NetCDF is a self-documenting binary format, making it easy to read and manipulate the output of a given module before feeding it to the next module. Furthermore, it is cross-platform meaning a processing chain can even be split to multiple different platforms, e.g. generating an L1A product using the instrument model on a Windows machine and using it as input to the L1A-L1B processor on a Linux machine.


Geometry module
------------------

The geometry module (GM) has a single driver function which takes a dictionary of user-defined parameters as input. Assuming TEDS is in your ``PYTHONPATH``, the following statements are sufficient to run the GM:

.. code-block:: python

   import yaml
   from teds.gm import geometry_module
   conf = yaml.safe_load(open('gm.yaml'))
   geometry_module(conf)

Here a GM configuration file is assumed to be in the working directory.


Scene generation module
-------------------------

The scene generation module (SGM) has three entry points:

.. code-block:: python

   import yaml
   from teds.sgm import Carbon_radiation_scene_generation
   from teds.sgm import geoscene_generation
   from teds.sgm import download_sentinel2_albedo
   conf = yaml.safe_load(open('sgm.yaml'))
   download_sentinel2_albedo(conf)
   geoscene_generation(conf)
   Carbon_radiation_scene_generation(conf)

The first step is to download a set of Sentinel 2 albedos and save them to a NetCDF file specified by ``[io_files][input_s2]`` in the configuration file. Unless the scene changes, this typically only needs to be done once. Next step is to generate the geophysical scene which defines the atmosphere. The final step is generating top of the atmosphere spectra for a given instrument. These serve as input to the instrument model.


Instrument model
------------------

The instrument model (IM) has both Python and C++ implementations. For the Python version, run

.. code-block:: python

   import yaml
   from teds.im import process_im
   conf = yaml.safe_load(open('im.yaml'))
   process_im(conf)

For the C++ version, if the build was successful, there is an executable called ``tango_im.x`` in the build directory. It takes one argument which is the YAML file that specifies how to apply the CKD to the line-by-line spectra produced by the scene generation module in order to generate a L1A product. It can be run with

.. code-block:: bash

   export OMP_NUM_THREADS=8
   <tango_im.x> im.yaml

The IM is parallelized over ALT positions (detector images) using OpenMP. If you exclude the ``export`` statement then the default is to run using all available threads.


L1A-L1B processor
------------------

The L1A-L1B processor is analogous to the IM in that it has both Python and C++ implementations. For the Python version, run

.. code-block:: python

   import yaml
   from teds.l1al1b import process_l1b
   conf = yaml.safe_load(open('l1b.yaml'))
   process_l1b(conf)

For the C++ version, run

.. code-block:: bash

   export OMP_NUM_THREADS=8
   <tango_l1b.x> l1b.yaml

The result is a level 1B product that can be used as input to the L1-L2 processor.


L1-L2 processor
----------------

A minimal script to run the L1-L2 processor is

.. code-block:: python

   import yaml
   from teds.l1l2 import level1b_to_level2_processor
   conf = yaml.safe_load(open('l2.yaml'))
   level1b_to_level2_processor(conf)


Full chain
----------

Having TEDS written as Python modules with no main driver function leaves the user with maximum flexibility in how they want to run the end-to-end chain. The following is a minimal example of how to run a full end-to-end simulation:

.. code-block:: bash

   # Import all relevant modules
   from subprocess import run  # For C++ modules
   import yaml
   from teds.ckd import gen_ckd
   from teds.gm import geometry_module
   from teds.sgm import Carbon_radiation_scene_generation
   from teds.sgm import geoscene_generation
   from teds.sgm import download_sentinel2_albedo
   from teds.im import process_im
   from teds.l1al1b import process_l1b
   from teds.siml1b import simplified_instrument_model_and_l1b_processor
   from teds.l1l2 import level1b_to_level2_processor

   # Read in all configuration files even if not required
   ckd_conf = yaml.safe_load(open('ckd.yaml'))
   gm_conf = yaml.safe_load(open('gm.yaml'))
   sgm_conf = yaml.safe_load(open('sgm.yaml'))
   siml1b_conf = yaml.safe_load(open('siml1b.yaml'))
   im_conf = yaml.safe_load(open('im.yaml'))
   l1b_conf = yaml.safe_load(open('l1b.yaml'))
   l2_conf = yaml.safe_load(open('l2.yaml'))

   # Run all or selected TEDS modules
   gen_ckd(ckd_conf)
   geometry_module(gm_conf)
   download_sentinel2_albedo(sgm_conf)
   geoscene_generation(sgm_conf)
   Carbon_radiation_scene_generation(sgm_conf)
   # Python IM & L1B
   process_im(im_conf)
   process_l1b(l1b_conf)
   # C++ IM & L1B
   # run(['tango_im.x', 'im.yaml'])
   # run(['tango_l1b.x', 'l1b.yaml'])
   level1b_to_level2_processor(l2_conf)

It is easy to comment out individual steps (as long as the chain remains consistent) or add control flow statements (loops, conditionals, ..) depending on the nature of the study (e.g. a sensitivity analysis). Note that here the statements are ordered such that the static parts (imports and configuration files) come first. Those can always be executed even if not needed. This leaves the part where the modules are run cleaner and easier to work with. Also, only the Python versions of the IM and L1A-L1B are run in this example. In order to run the C++ versions, uncomment them. The C++ executables are assumed to be in ``PATH`` here. If not, provide the full paths instead.
