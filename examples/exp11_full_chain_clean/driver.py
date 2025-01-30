# Example of Tango Carbon full-chain simulation. This is a "clean"
# version, meaning no complex path or logic dependencies. The driver
# script sequentially runs all TEDS modules by reading and writing to
# the present directory. The configuration files only contain relative
# paths assuming all files are in the present directory. The auxiliary
# files can be downloaded from
# https://surfdrive.surf.nl/files/index.php/s/xs0PS5KyxvNxYTX. For
# this example, the following files and directories are required:
#   hapi
#   binning_table.nc
#   ckd_stray_1.nc
#   co2_src1_20180523_1230.nc
#   gebco_ocssw_v2020.nc
#   OWL640S_low-gain_ckd.nc
#   prof.AFGL.US.std
#   solar_spectrum.nc
#   spot_col_distances.dat
#   spot_row_distances.dat
#
# Run as follows:
#
#  # Activate virtual environment
#  # Download auxiliary files into the present directory
#  python driver.py
#
# This example excludes the simplified L1B processor and runs the
# Python IM and L1B. The C++ IM and L1B are included in the list of
# modules but commented out below. In order to run them, the C++
# modules must be compiled first and the executables should either be
# in PATH or the executable paths need to be adjusted accordingly
# below.

# from subprocess import run  # For running the C++ IM and L1B
import yaml

from teds.ckd import gen_ckd
from teds.gm.gm import geometry_module
from teds.sgm import geoscene_generation
from teds.sgm import Carbon_radiation_scene_generation
from teds.sgm.s2 import download_albedo
from teds.im import run_instrument_model
from teds.l1al1b import run_l1al1b
from teds.l1l2 import level1b_to_level2_processor

gen_ckd(yaml.safe_load(open('ckd.yaml')))
geometry_module(yaml.safe_load(open('gm.yaml')))
download_albedo(yaml.safe_load(open('sgm.yaml')))
geoscene_generation(yaml.safe_load(open('sgm.yaml')))
Carbon_radiation_scene_generation(yaml.safe_load(open('sgm.yaml')))
run_instrument_model(yaml.safe_load(open('im.yaml')))
run_l1al1b(yaml.safe_load(open('l1b.yaml')))
# run(['tango_im.x', 'im.yaml'])  # Assuming the executable is in PATH
# run(['tango_l1b.x', 'l1b.yaml'])  # Assuming the executable is in PATH
level1b_to_level2_processor(yaml.safe_load(open('l2.yaml')))
