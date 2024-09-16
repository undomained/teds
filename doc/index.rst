============================================
Toolbox for End-to-enD Simulations (TEDS)
============================================

v |version|

TEDS is a toolbox of harmonized software modules for end-to-end (E2E) performance analyses for the Twin ANthropogenic Greenhouse gas Observers (TANGO) mission. TANGO is a SCOUT mission within the ESA's Earth Observation FutureEO program to quantify and monitor greenhouse gas emissions at the level of individual facilities using two CubeSat satellites.

The TEDS package includes:
 * Geometry module -- satellite ephemeris and attitude
 * Scene generation module -- simulate atmosphere and a VIS and SWIR radiation scene
 * Instrument model -- convert radiation scene into raw detector images (level 0)
 * L1A-L1B processor -- convert detector images into radiation scene (level 1B)
 * L1B-L2 processor -- convert L1B product to CO\ :sub:`2`\ , CH\ :sub:`4`\ , and NO\ :sub:`2` columns (level 2)
 * L2-L4 processor -- plume simulation (level 4)

The TEDS software modules are designed to be used in an E2E software pipeline with file-based interfaces between the modules. For complete mission performance analyses, platform and instrument specifications need to be provided as input to TEDS modules.

This documentation explains how to build or install and run different modules of TEDS. The algorithms are described in correspsonding Algorithm Theoretical Baseline Documents (ATBDs).


Contents
---------

.. toctree::
   :maxdepth: 2

   downloading
   installation
   running_the_code
   api
   contributing
