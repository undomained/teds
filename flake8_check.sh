#!/bin/bash

# Run flake8 to check for PEP 8 violations on Python source files
# including tests. Some source files are exluded from the check for
# now. No new exceptions should be added.

files=$(find teds -name "*.py" -regextype posix-extended \
             ! -wholename teds/lib/hapi.py \
             ! -wholename teds/ckd/ckd_generation/ckd_generator_nitro.py \
             ! -wholename teds/ckd/ckd_generation/generator_class.py \
             ! -wholename teds/ckd/create_binning_tables.py \
             ! -wholename "teds/ckd/ckd_generation/nitro_ckds/*" \
             ! -wholename teds/sgm/sgm_no2.py \
             ! -wholename teds/im/create_im_configuration_file_nitro.py \
             ! -wholename teds/l1al1b/create_l1a1b_configuration_file_nitro.py \
             ! -wholename teds/l1l2/l1bl2_no2.py \
             ! -wholename teds/lib/libRT_no2.py \
             ! -wholename teds/lib/libDOAS.py \
             ! -wholename teds/lib/libAMF.py \
             ! -wholename "teds/lib/data_netcdf/*" \
         )
flake8 $files
exit $?
