# Find oasis
#
# Please pass the -DOASIS_DIR variable to cmake (location of the OASIS libraries).

# Search for all library files individually and merge to a single OASIS_LIBRARIES variable for the parent CMakeLists.txt
find_library(LIB_OASIS NAMES oasis.cbind HINTS ${OASIS_DIR}/lib ENV LD_LIBRARY_PATH)
find_library(LIB_MCT NAMES mct HINTS ${OASIS_DIR}/lib ENV LD_LIBRARY_PATH)
find_library(LIB_MPEU NAMES mpeu HINTS ${OASIS_DIR}/lib ENV LD_LIBRARY_PATH)
find_library(LIB_PSMILE NAMES psmile.MPI1 HINTS ${OASIS_DIR}/lib ENV LD_LIBRARY_PATH)
find_library(LIB_SCRIP NAMES scrip HINTS ${OASIS_DIR}/lib ENV LD_LIBRARY_PATH)

if (LIB_OASIS AND LIB_MCT AND LIB_MPEU AND LIB_PSMILE AND LIB_SCRIP)
    set(OASIS_LIBRARIES ${LIB_OASIS} ${LIB_MCT} ${LIB_MPEU} ${LIB_PSMILE} ${LIB_SCRIP})
endif ()

get_filename_component(OASIS_LIB_DIR "${LIB_OASIS}" PATH)

set(OASIS_DIR "${OASIS_LIB_DIR}/../")
cmake_path(NORMAL_PATH OASIS_DIR)

find_path(OASIS_INCLUDES NAMES oasis_c.h HINTS ${OASIS_DIR}/include)

# handle the QUIETLY and REQUIRED arguments and set OASIS_FOUND to TRUE if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OASIS DEFAULT_MSG OASIS_DIR OASIS_LIBRARIES OASIS_INCLUDES)

# message (STATUS "OASIS_LIBRARIES: ${OASIS_LIBRARIES}")
# message (STATUS "OASIS_INCLUDES: ${OASIS_INCLUDES}")
# message (STATUS "OASIS_DIR: ${OASIS_DIR}")
