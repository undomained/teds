#!/bin/bash

# Run flake8 to check for PEP 8 violations on all Python source files,
# including the tests.

all_files=$(find teds -name "*.py" -regextype posix-extended -not -regex "teds/(lib/hapi|lib/data_netcdf/.*|lib/libAMF|lib/libDOAS|lib/libRT_no2|sgm/sgm_no2|im/create_im_configuration_file_nitro|ckd/ckd_generation/ckd_generator_nitro|ckd/ckd_generation/generator_class|ckd/ckd_generation/nitro_ckds/.*|ckd/create_binning_tables|l1al1b/create_l1a1b_configuration_file_nitro|l1l2/l1bl2_no2|siml1b/radiance_offset).py")
flake8 $all_files
exit $?
