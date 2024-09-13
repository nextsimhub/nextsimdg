# Find oasis
#
# Please pass the -DOASIS_DIR variable to cmake (location of the OASIS libraries).

find_library (OASIS_LIBRARIES NAMES mct mpeu oasis psmile scrip HINTS ${OASIS_DIR}/lib ENV LD_LIBRARY_PATH)

get_filename_component (OASIS_LIBRARIES "${OASIS_LIBRARIES}" PATH)

set (OASIS_DIR "${OASIS_LIBRARIES}/../")
cmake_path(NORMAL_PATH OASIS_DIR)

find_path (OASIS_INCLUDES NAMES oasis_c.h HINTS ${OASIS_DIR}/include)

# handle the QUIETLY and REQUIRED arguments and set OASIS_FOUND to TRUE if all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(OASIS DEFAULT_MSG
        OASIS_LIBRARIES OASIS_INCLUDES)

# message (STATUS "OASIS_LIBRARIES: ${OASIS_LIBRARIES}")
# message (STATUS "OASIS_INCLUDES: ${OASIS_INCLUDES}")
# message (STATUS "OASIS_DIR: ${OASIS_DIR}")
