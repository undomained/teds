# This source code is licensed under the 3-clause BSD license found in
# the LICENSE file in the root directory of this project.

if (__GET_GIT_INFO_GUARD__)
  return()
endif()
set(__GET_GIT_INFO_GUARD__ TRUE)

# Get full hash of the latest git commit
execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND git log --format=%H -1
  OUTPUT_VARIABLE GIT_COMMIT)
string(STRIP ${GIT_COMMIT} GIT_COMMIT)
# Find out if there are uncommited changes
execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND git status --untracked-files=no --porcelain
  OUTPUT_VARIABLE GIT_UNTRACKED)
if (GIT_UNTRACKED)
  set(GIT_COMMIT "${GIT_COMMIT} with changes")
endif()
# Get date of the latest commit
execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  COMMAND git log --format=%aD -1
  OUTPUT_VARIABLE GIT_DATE)
string(REGEX REPLACE "(^[^ ]| [^ ]+ [^ ]+$)" "" GIT_DATE ${GIT_DATE})
string(STRIP ${GIT_DATE} GIT_DATE)
# Get the latest tag
# execute_process(WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
#   COMMAND git describe --tags
#   OUTPUT_VARIABLE GIT_TAG)
#string(REGEX REPLACE "-.+$" "" GIT_TAG ${GIT_TAG})
#string(REGEX REPLACE "\n" "" GIT_TAG ${GIT_TAG})
configure_file(${PROJECT_SOURCE_DIR}/tango/git.h.in git.h)
