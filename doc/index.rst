=========================================
Tango end-to-end simulator v |version|
=========================================

Tango end-to-end simulator (TEDS) provides a set of tools to simulate level 0-4 data products for the Twin ANthropogenic Greenhouse gas Observers (Tango) mission. Tango is a SCOUT mission within the ESA's Earth Observation FutureEO program. It will quantify and monitor greenhouse gas emissions at the level of individual facilities using two CubeSat satellites.

This documentation explains how to build or install and run different modules of the TEDS package:
 * Geometry module -- satellite ephemeris and attitude
 * Scene generation module -- simulate atmosphere and a VIS and SWIR radiation scene
 * Instrument model -- convert radiation scene into raw detector images
 * L1A-L1B processor -- convert detector images into radiation scene
 * L1B-L2 processor -- convert L1B product to CO\ :sub:`2`\ , CH\ :sub:`4`\ , and NO\ :sub:`2` columns (L2 product)
 * L2-L4 processor -- plume simulation

Contents
---------

.. toctree::
   :maxdepth: 2

   downloading
   installation
   running_the_code
   contributing

Terms of use
-------------

TEDS is licensed under the 3-clause BSD license found in the LICENSE file in the root directory of the project.

Contributors
-------------

In order to see a list of contributors, clone the repository and issue

.. code-block:: bash

   git shortlog -s -n --all --no-merges

This prints a list of contributors along with a number of commits by each contributor.
