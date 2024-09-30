===
API
===

This section provides an API for Python-based modules. Each section lists the main entry points first, followed by functions describing the internals of each module (mostly for developers).


Scene generation module (SGM)
-----------------------------

.. autofunction:: teds.sgm.download_sentinel2_albedo
.. autofunction:: teds.sgm.geoscene_generation
.. autofunction:: teds.sgm.Carbon_radiation_scene_generation


Instrument model
-----------------

Main entry point
^^^^^^^^^^^^^^^^

.. autofunction:: teds.im.process_im

Uncalibration steps
^^^^^^^^^^^^^^^^^^^

.. autofunction:: teds.im.python.uncalibration.change_wavelength_grid
.. autofunction:: teds.im.python.uncalibration.convert_from_radiance
.. autofunction:: teds.im.python.uncalibration.map_to_detector
.. autofunction:: teds.im.python.uncalibration.stray_light
.. autofunction:: teds.im.python.uncalibration.include_prnu
.. autofunction:: teds.im.python.uncalibration.include_nonlinearity
.. autofunction:: teds.im.python.uncalibration.include_dark_signal
.. autofunction:: teds.im.python.uncalibration.include_noise
.. autofunction:: teds.im.python.uncalibration.include_offset
.. autofunction:: teds.im.python.uncalibration.include_coadding_and_binning


Input/output
^^^^^^^^^^^^

These are the same as for the L1B processor below.


L1B processor
-------------

Main entry point
^^^^^^^^^^^^^^^^

.. autofunction:: teds.l1al1b.process_l1b


Calibration steps
^^^^^^^^^^^^^^^^^

.. autofunction:: teds.l1al1b.python.calibration.remove_coadding_and_binning
.. autofunction:: teds.l1al1b.python.calibration.remove_offset
.. autofunction:: teds.l1al1b.python.calibration.determine_noise
.. autofunction:: teds.l1al1b.python.calibration.remove_dark_signal
.. autofunction:: teds.l1al1b.python.calibration.remove_nonlinearity
.. autofunction:: teds.l1al1b.python.calibration.remove_prnu
.. autofunction:: teds.l1al1b.python.calibration.stray_light
.. autofunction:: teds.l1al1b.python.calibration.map_from_detector
.. autofunction:: teds.l1al1b.python.calibration.convert_to_radiance

Input/output
^^^^^^^^^^^^

.. autofunction:: teds.l1al1b.python.io.read_proc_level
.. autofunction:: teds.l1al1b.python.io.read_ckd
.. autofunction:: teds.l1al1b.python.io.read_geometry
.. autofunction:: teds.l1al1b.python.io.read_l1
