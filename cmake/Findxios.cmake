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

set(xios_INCLUDES "${xios_DIR}/inc/")
set(xios_LIBRARIES "${xios_DIR}/lib/")
set(xios_EXTERNS "${xios_DIR}/extern/")

find_path (xios_INCLUDES xios.hpp
  HINTS xios_DIR ENV xios_DIR ENV CPATH ENV FPATH)

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

find_path (xios_INCLUDES xios.hpp
  HINTS xios_DIR ENV xios_DIR ENV CPATH ENV FPATH)

find_library (xios_LIBRARIES NAMES xios
    HINTS ENV LD_LIBRARY_PATH LIBRARY_PATH)
mark_as_advanced(xios_LIBRARIES)

set (xois_has_interfaces "YES") # will be set to NO if we're missing any interfaces
set (xios_libs "${xios_LIBRARIES}")

get_filename_component (xios_lib_dirs "${xios_LIBRARIES}" PATH)

set (xios_LIBRARIES "${xios_libs}" CACHE STRING "All xios libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set xios_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(xios DEFAULT_MSG xios_LIBRARIES xios_INCLUDES)# xios_has_interfaces)

mark_as_advanced(xios_LIBRARIES xios_INCLUDES xios_EXTERNS)

message (STATUS "xios_LIBRARIES: ${xios_LIBRARIES}")
message (STATUS "xios_INCLUDES: ${xios_INCLUDES}")
