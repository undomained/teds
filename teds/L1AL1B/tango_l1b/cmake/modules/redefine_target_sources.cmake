if (__REDEFINE_TARGET_SOURCES_GUARD__)
  return()
endif()
set(__REDEFINE_TARGET_SOURCES_GUARD__ TRUE)

# Emulate the current behavior of target_sources for older CMake
# versions.
macro(target_sources target_name mode)
  foreach(_src ${ARGN})
    if (EXISTS ${_src})
      _target_sources(${target_name} ${mode} ${_src})
    else()
      _target_sources(${target_name} ${mode} ${CMAKE_CURRENT_LIST_DIR}/${_src})
    endif()
  endforeach()
endmacro()
