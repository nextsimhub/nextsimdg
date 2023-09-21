# Find xios
# 
# 26th January 2023
# Dr Alexander Smith, ICCS, Cambridge (UK)
#
# Find the xios includes and library 
#
# Please set the xios_DIR to point to the top level directory for xios
# (which contains the makexios executable). This is <path>/xios/trunk/
# if you followed the instructions in the user/reference guide.
#

#if (NOT (xios_DIR))
  set( xios_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xios/")
set(xios_INCLUDES "${xios_DIR}/inc/")
set(xios_LIBRARIES "${xios_DIR}/lib/")
set(xios_EXTERNS "${xios_DIR}/extern/")

find_path (xios_INCLUDES xios.hpp
  HINTS xios_DIR ENV xios_DIR ENV CPATH ENV FPATH)
#set(CMAKE_FIND_DEBUG_MODE FALSE)
message (STATUS "so far xios_INCLUDES: ${xios_INCLUDES}")

# - Find NetCDF
# Find the native NetCDF includes and library
#
#  xios_INCLUDES    - where to find netcdf.h, etc
#  xios_LIBRARIES   - Link these libraries when using NetCDF
#  xios_FOUND       - True if NetCDF found including required interfaces (see below)

#  find_package (xios REQUIRED)
#  target_link_libraries (uses_f90_interface ${xios_LIBRARIES})
#  target_link_libraries (only_uses_c_interface ${xios_LIBRARIES_C})

if (xios_INCLUDES AND xios_LIBRARIES)
  # Already in cache, be silent
  set (xios_FIND_QUIETLY TRUE)
endif (xios_INCLUDES AND xios_LIBRARIES)

#set(CMAKE_FIND_DEBUG_MODE TRUE)
find_path (xios_INCLUDES xios.hpp
  HINTS xios_DIR ENV xios_DIR ENV CPATH ENV FPATH)
#set(CMAKE_FIND_DEBUG_MODE FALSE)
message (STATUS "so far xios_INCLUDES: ${xios_INCLUDES}")


find_library (xios_LIBRARIES NAMES xios
    HINTS ENV LD_LIBRARY_PATH LIBRARY_PATH)
mark_as_advanced(xios_LIBRARIES)
message (STATUS "so far xios_LIBRARIES: ${xios_LIBRARIES}")

set (xois_has_interfaces "YES") # will be set to NO if we're missing any interfaces
set (xios_libs "${xios_LIBRARIES}")

get_filename_component (xios_lib_dirs "${xios_LIBRARIES}" PATH)

message (STATUS "lang is set to:  ${lang}")

#macro (NetCDF_check_interface lang header libs)
#  if (xios_${lang})
#    find_path (xios_INCLUDES_${lang} NAMES ${header}
#      HINTS "${xios_INCLUDES}" NO_DEFAULT_PATH)
#    find_library (xios_LIBRARIES_${lang} NAMES ${libs}
#      HINTS "${NetCDF_lib_dirs}" NO_DEFAULT_PATH)
#    mark_as_advanced (xios_INCLUDES_${lang} xios_LIBRARIES_${lang})
#    if (xios_INCLUDES_${lang} AND xios_LIBRARIES_${lang})
#      list (INSERT NetCDF_libs 0 ${xios_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
#    else (xios_INCLUDES_${lang} AND xios_LIBRARIES_${lang})
#      set (NetCDF_has_interfaces "NO")
#      message (STATUS "Failed to find NetCDF interface for ${lang}")
#    endif (xios_INCLUDES_${lang} AND xios_LIBRARIES_${lang})
#  endif (xios_${lang})
#endmacro (NetCDF_check_interface)

#NetCDF_check_interface (CXX netcdfcpp.h netcdf_c++)
#NetCDF_check_interface (F77 netcdf.inc  netcdff)
#NetCDF_check_interface (F90 netcdf.mod  netcdff)

set (xios_LIBRARIES "${xios_libs}" CACHE STRING "All xios libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set xios_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(xios DEFAULT_MSG xios_LIBRARIES xios_INCLUDES)# xios_has_interfaces)

mark_as_advanced(xios_LIBRARIES xios_INCLUDES xios_EXTERNS)


