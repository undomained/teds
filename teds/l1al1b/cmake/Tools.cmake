if (__TOOLS_GUARD__)
  return()
endif()
set(__TOOLS_GUARD__ TRUE)

# Convert a string into a list. If the variable was already given as
# list, do nothing.
macro(convert_to_list var)
  list(LENGTH ${var} _length)
  if (_length EQUAL 1)
    string(REPLACE " " ";" ${var} ${${var}})
  endif()
endmacro()

# Go through directories listed in LIB_PATHS and turn the entries of
# LIBS into actual targets.
function(generate_library_targets _PATHS _LIBS)
  foreach(LIB ${${_LIBS}})
    # If the target already exists, skip it. This can occur when the
    # same library is set multiple times in LIBS due to circular
    # dependencies.
    if (NOT TARGET ${LIB})
      find_library(LIB_FULLPATH ${LIB} HINTS ${${_PATHS}})
      if (LIB_FULLPATH)
        message(STATUS "Found ${LIB_FULLPATH}")
        add_library(${LIB} UNKNOWN IMPORTED)
        set_target_properties(${LIB} PROPERTIES
          IMPORTED_LOCATION ${LIB_FULLPATH})
        unset(LIB_FULLPATH CACHE)
      else()
        message(FATAL_ERROR "Could not find ${LIB}")
      endif()
    endif()
  endforeach()
endfunction()

# Determine current flags
macro(get_all_flags language flags)
  set(${flags} ${CMAKE_${language}_FLAGS})
  if (CMAKE_BUILD_TYPE)
    string(TOUPPER ${CMAKE_BUILD_TYPE} buildtype)
    set(${flags} "${CMAKE_${language}_FLAGS_${buildtype}} ${${flags}}")
    string(STRIP ${${flags}} ${flags})
  endif()
endmacro()

# Return a list of paths of libraries linked to a target. The return
# value (${target}_library_paths) is a list of strings.
function(get_library_paths target)
  get_target_property(_libraries ${target} LINK_LIBRARIES)
  set(_paths "")
  foreach(_lib ${_libraries})
    # Test if any of the following properties are defined. Those are
    # all candidates for the where the library is located.
    get_target_property(_loc ${_lib} IMPORTED_LOCATION)
    if (_loc MATCHES ".*-NOTFOUND")
      get_target_property(_loc ${_lib} IMPORTED_LOCATION_RELEASE)
    endif()
    if (_loc MATCHES ".*-NOTFOUND")
      get_target_property(_loc ${_lib} INTERFACE_INCLUDE_DIRECTORIES)
    endif()
    # If non of the above was found simply set the location to the
    # library name (basically remove the NOTFOUND part).
    if (_loc MATCHES ".*-NOTFOUND")
      set(_loc ${_lib})
    endif()
    # Clean the path if from a build instead of an install directory
    string(REGEX REPLACE "\\$<BUILD_INTERFACE:([^>]+).+" "\\1" _loc "${_loc}")
    list(APPEND _paths ${_loc})
  endforeach()
  list(JOIN _paths " " _paths)
  set(${target}_library_paths ${_paths} PARENT_SCOPE)
endfunction()

# Test if INCLUDE_COVERAGE option is ON and if so pass the coverage
# flag to the compiler. This macro needs to be called separately for
# each scope.
macro(check_include_coverage)
  if (INCLUDE_COVERAGE)
    string(APPEND CMAKE_CXX_FLAGS " ${COVERAGE_FLAGS}")
  endif()
endmacro()
