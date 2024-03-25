# Find xios
#
# Please set the xios_DIR to point to the top level directory for xios
# (which contains the makexios executable). This is <path>/xios/trunk/
# if you followed the instructions in the user/reference guide.
#

find_library (xios_LIBRARIES NAMES xios HINTS ${xios_DIR}/lib ENV LD_LIBRARY_PATH)

get_filename_component (xios_LIBRARIES "${xios_LIBRARIES}" PATH)

set (xios_DIR "${xios_LIBRARIES}/../")
cmake_path(NORMAL_PATH xios_DIR)

find_path (xios_INCLUDES NAMES xios.hpp HINTS ${xios_DIR}/inc)

set (xios_EXTERNS "${xios_DIR}/extern")
cmake_path(NORMAL_PATH xios_EXTERNS)

# handle the QUIETLY and REQUIRED arguments and set xios_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(xios DEFAULT_MSG
        xios_LIBRARIES xios_INCLUDES xios_EXTERNS)

message (STATUS "xios_LIBRARIES: ${xios_LIBRARIES}")
message (STATUS "xios_INCLUDES: ${xios_INCLUDES}")
message (STATUS "xios_EXTERNS: ${xios_EXTERNS}")
# message (STATUS "xios_DIR: ${xios_DIR}")
