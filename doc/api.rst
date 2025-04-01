===
API
===

This section provides an API for Python-based modules. Each section lists the main entry points first, followed by functions describing the internals of each module (mostly for developers).


Scene generation module (SGM)
-----------------------------

.. autofunction:: teds.sgm.s2.download_albedo
.. autofunction:: teds.sgm.geoscene_generation
.. autofunction:: teds.sgm.carbon_radiation_scene_generation


Instrument model
-----------------

Main entry point
^^^^^^^^^^^^^^^^

.. autofunction:: teds.im.run_instrument_model

Forward model
^^^^^^^^^^^^^^

.. autofunction:: teds.im.python.forward_models.apply_isrf
.. autofunction:: teds.im.python.forward_models.radiometric
.. autofunction:: teds.im.python.forward_models.map_to_detector
.. autofunction:: teds.im.python.forward_models.stray_light
.. autofunction:: teds.im.python.forward_models.prnu
.. autofunction:: teds.im.python.forward_models.nonlinearity
.. autofunction:: teds.im.python.forward_models.dark_current
.. autofunction:: teds.im.python.forward_models.noise
.. autofunction:: teds.im.python.forward_models.dark_offset
.. autofunction:: teds.im.python.forward_models.coadd_and_adc


Input/output
^^^^^^^^^^^^

These are the same as for the L1B processor below.


L1B processor
-------------

Main entry point
^^^^^^^^^^^^^^^^

.. autofunction:: teds.l1al1b.run_l1al1b


Calibration steps
^^^^^^^^^^^^^^^^^

.. autofunction:: teds.l1al1b.python.calibration.coadding_and_binning
.. autofunction:: teds.l1al1b.python.calibration.dark_offset
.. autofunction:: teds.l1al1b.python.calibration.noise
.. autofunction:: teds.l1al1b.python.calibration.dark_current
.. autofunction:: teds.l1al1b.python.calibration.nonlinearity
.. autofunction:: teds.l1al1b.python.calibration.prnu
.. autofunction:: teds.l1al1b.python.calibration.stray_light
.. autofunction:: teds.l1al1b.python.calibration.map_from_detector
.. autofunction:: teds.l1al1b.python.calibration.radiometric

Input/output
^^^^^^^^^^^^

.. autofunction:: teds.l1al1b.python.io.read_proc_level
.. autofunction:: teds.l1al1b.python.io.read_ckd
.. autofunction:: teds.l1al1b.python.io.read_geometry
.. autofunction:: teds.l1al1b.python.io.read_l1
