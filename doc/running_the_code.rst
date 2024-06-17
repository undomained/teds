Running the code
====================

If the build process was successful, an executable called ``tango_l1b.x`` is produced in the build directory ``<build>``. The executable takes one argument which is a configuration file that specifies how to apply the CKD to flight images:

.. code-block:: bash

   <build>/tango_l1b.x run.yaml

where ``run.yaml`` is a user defined configuration file. If the processor was compiled using a parallel C++ compiler, it can be run by issuing

.. code-block:: bash

   export OMP_NUM_THREADS=8
   <build>/tango_l1b.x run.yaml


YAML configuration files
-------------------------

Tango configuration files are written in YAML. We chose YAML because it is versatile, human-readable, and widely used and supported. It uses Python-style indentation to indicate nesting and as basic types supports scalars, lists, and maps (also called dictionaries or hashes). Furthermore, it is a strict superset of JSON and thus compatible with it.

In order to understand the YAML format open one of the tutorial configuration files or run the processor without any arguments:

.. code-block:: bash

   <build>/tango_l1b.x > conf.yaml

The output contains keywords showing their type, default value, and description. Many default values are not suitable and must be edited before feeding the configuration file to the processor. The configuration file may contain several sections. For example, a ``detector`` section contains a set of parameters specifically targeting the detector.

The basic element of YAML is a key-value pair where the value is a scalar, list, or a map. Here is an example where the element type is noted as a comment::

  main: # starts a map
    processes: # key followed by a list
      - pol # scalar
      - l1b # scalar
  ckd_file_in: ckd.nc # key-value
  detector_spec: 2048 # key-value

When you're done editing the configuration file, run the processor with the argument:

.. code-block:: bash

   <build>/tango_l1b.x conf.yaml

There are many checks in the code which notify you of errors and inconsistencies in the configuration file. If you use a non-existent keyword the processor will issue a warning but will continue with the calculation.

Full specification of YAML and how to use it is found here: https://yaml.org/spec/1.2.2/.
