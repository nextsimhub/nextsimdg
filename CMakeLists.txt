cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)

project(nextsim_dg CXX)

# If OASIS is requested, then Fortran needs to be enabled for the project.
# Take this opportunity to define the basic OASIS variables.
if(DEFINED OASIS_LIB_DIR)
    add_compile_definitions(USE_OASIS)
    set(USE_OASIS TRUE)
    enable_language(Fortran)
else()
    set(USE_OASIS FALSE)
endif()

find_package(PkgConfig)

# OpenMP and MPI libraries
find_package(OpenMP REQUIRED)
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

find_package(MPI REQUIRED)
set(NSDG_MPI_Libs MPI::MPI_CXX)

pkg_search_module(NETCDF_CXX4 netcdf-cxx4)
if (NETCDF_CXX4_FOUND)
    set(NSDG_NetCDF_Library "${NETCDF_CXX4_LIBRARIES}")
    set(netCDF_INC_DIR "${NETCDF_CXX4_INCLUDE_DIRS}")
    set(netCDF_LIB_DIR "${NETCDF_CXX4_LIBRARY_DIRS}")
else()
    find_package(netCDF REQUIRED)
    if("${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Darwin")
        set(NSDG_NetCDF_Library "netcdf-cxx4")
    else()
        set(NSDG_NetCDF_Library "netcdf_c++4")
    endif()
endif()
# To add netCDF to a target:
# target_include_directories(target PUBLIC ${netCDF_INCLUDE_DIR})
# target_link_directories(target PUBLIC ${netCDF_LIB_DIR})
# target_link_libraries(target LINK_PUBLIC "${NSDG_NetCDF_Library}")

# Regarding Boost.Log, if our application consists of more
# than one modules that use it, we must link to the shared
# version. If we have a single executable or a single module
# that works, we may use the static version.
# By default, it is assumed that the library is built in
# static mode. Use the following definition to indicate that
# the code will be linked against dynamically loaded boost
# libraries.
add_definitions(-DBOOST_ALL_DYN_LINK)
find_package(Boost COMPONENTS program_options log REQUIRED)

find_package(Catch2 REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

# OASIS libraries. The OASIS_LIB_DIR should be passed to activate OASIS linking
if(USE_OASIS)
    set(CHAN "MPI1")
    set(OasisLibraryDir "${OASIS_LIB_DIR}")
    set(PsmileLibraries
        "libpsmile.${CHAN}.a"
        "libmct.a"
        "libmpeu.a"
        "libscrip.a"
        "liboasis.cbind.a"
    )
    set(OasisIncludeDir "${OasisLibraryDir}/cbindings/c_src")
    set(NSDG_MPI_Libraries "${NSDG_MPI_Libs}" MPI::MPI_Fortran)
else()
    # Define an empty library to keep all the OASIS conditionals in one place
    set(PsmileLibraries "")
    set(OasisIncludeDir "")
endif()

# Set the list of components that define modules. Component subdirectories are
# contained in this directory level and contain src/modules/ and 
# src/modules/include/ subdirectories 
set(ModularComponents
"physics"
)

# Set the components which provide source or object code to the main model
set(CodeComponents
"physics"
#"dynamics"
)

# Set an empty list of sources
set(NextsimSources "")

# Set an empty list of include directories
set(NextsimIncludeDirs "")

# Build the core model. Defines the 'parse_modules' target
add_subdirectory(core)

# Build all components
foreach(compo ${CodeComponents})
    add_subdirectory("${compo}")
endforeach()

# Link the model
add_executable(nextsim "${NextsimSources}")

target_include_directories(nextsim PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NextsimIncludeDirs}"
    "${netCDF_INC_DIR}"
    "${OasisIncludeDir}"
    )
target_link_directories(nextsim PUBLIC "${netCDF_LIB_DIR}" "${OasisLibraryDir}")
target_link_libraries(nextsim LINK_PUBLIC
    Boost::program_options
    Boost::log
    "${NSDG_NetCDF_Library}"
    Eigen3::Eigen
    "${NSDG_MPI_Libraries}"
    "${PsmileLibraries}"
    )
