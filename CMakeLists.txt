cmake_minimum_required(VERSION 3.12)

# set globally because we have a number of independent targets
set(CMAKE_CXX_STANDARD 17)

project(nextsim_dg)

# user definable build options
option(WITH_THREADS "Build with support for openmp" OFF)
option(ENABLE_MPI "Enable distributed-memory parallelization with MPI" OFF)
option(ENABLE_XIOS "Enable XIOS library for IO" OFF)
option(ENABLE_OASIS "Enable the OASIS interface" OFF)
option(BUILD_TESTS "Build the tests" ON)
set(DynamicsType
    "DG2"
    CACHE STRING
    "Set the discretization order for the dynamics, either DG1 or DG2"
)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# all nextsim code is collected in a library to allow for reuse in the tests
# the code is added after processing the subdirectories
add_library(nextsimlib SHARED)
# modern way to set the standard is redundant since we also do it globally
target_compile_features(nextsimlib PUBLIC cxx_std_17)

# NetCDF
find_package(PkgConfig)
pkg_search_module(NETCDF_CXX4 netcdf-cxx4)
if(NETCDF_CXX4_FOUND)
    set(NSDG_NetCDF_Library "${NETCDF_CXX4_LIBRARIES}")
    set(netCDF_INCLUDE_DIR "${NETCDF_CXX4_INCLUDE_DIRS}")
    set(netCDF_LIB_DIR "${NETCDF_CXX4_LIBRARY_DIRS}")
else()
    find_package(netCDF REQUIRED)
    if("${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Darwin")
        set(NSDG_NetCDF_Library "netcdf-cxx4")
    else()
        set(NSDG_NetCDF_Library "netcdf_c++4")
    endif()
endif()
target_include_directories(nextsimlib PUBLIC "${netCDF_INCLUDE_DIR}")
target_link_directories(nextsimlib PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(nextsimlib PUBLIC "${NSDG_NetCDF_Library}" "${FORTRAN_RUNTIME_LIB}")

if(BUILD_TESTS)
    # Add the doctest header library
    set(DOCTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
    add_library(doctest::doctest IMPORTED INTERFACE)
    set_property(
        TARGET doctest::doctest
        PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${DOCTEST_INCLUDE_DIR}"
    )
endif()

# Eigen
find_package(Eigen3 3.4 REQUIRED)
target_link_libraries(nextsimlib PUBLIC Eigen3::Eigen)

set(Python_FIND_VIRTUALENV FIRST)
find_package(Python COMPONENTS Interpreter)

# The location of the module_builder.py scripts
set(ScriptDirectory "${PROJECT_SOURCE_DIR}/scripts")
set(ModuleBuilderScript "${ScriptDirectory}/module_builder.py")
set(ModuleHeaderScript "${ScriptDirectory}/module_header.py")

# Dynamics type. Defaults to DG1.
# Available options are:
# DG1 Discontinuous Galerkin, degree 1 (3 components)
# DG2 Discontinuous Galerkin, degree 2 (6 components)
set(isDG FALSE)

message("Dynamics is set to ${DynamicsType}")

if(DynamicsType STREQUAL "DG1")
    set(isDG TRUE)
    set(DGComp 3)
    set(CGDegree 1)
endif()

if(DynamicsType STREQUAL "DG2")
    set(isDG TRUE)
    set(DGComp 6)
    set(CGDegree 2)
endif()

# Set the components which provide source or object code to the main model
set(CodeComponents "physics")

if(isDG)
    # Add the DG dynamics subdirectory to the list of code components
    set(CodeComponents "${CodeComponents}" "dynamics")

    # Set the number of DG stress components given the CG degree
    if(CGDegree EQUAL 1)
        set(DGStressComp 3)
    elseif(CGDegree EQUAL 2)
        set(DGStressComp 8)
    else()
        message("Invalid value of CGDegree. Valid values are 1–2")
    endif()

    # The exact selection of dimensions and Types available to ModelArray
    set(ModelArrayStructure "discontinuousgalerkin")
endif()
# set degree / component parameters
target_compile_definitions(
    nextsimlib
    PRIVATE DGCOMP=${DGComp} DGSTRESSCOMP=${DGStressComp} CGDEGREE=${CGDegree}
)

# boost
find_package(Boost COMPONENTS program_options log REQUIRED)
target_link_libraries(nextsimlib PUBLIC Boost::program_options Boost::log)
# Regarding Boost.Log, if our application consists of more
# than one modules that use it, we must link to the shared
# version. If we have a single executable or a single module
# that works, we may use the static version.
# By default, it is assumed that the library is built in
# static mode. Use the following definition to indicate that
# the code will be linked against dynamically loaded boost
# libraries.
target_compile_definitions(nextsimlib PUBLIC BOOST_ALL_DYN_LINK)

# MPI
if(ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS C CXX)
    target_compile_definitions(nextsimlib PUBLIC USE_MPI)
    target_link_libraries(nextsimlib PUBLIC MPI::MPI_C MPI::MPI_CXX)
endif()

# XIOS
if(ENABLE_XIOS)
    message(STATUS "Building with XIOS library")
    find_package(xios REQUIRED)
    find_package(HDF5 REQUIRED COMPONENTS C HL)
    target_compile_definitions(nextsimlib PUBLIC USE_XIOS)

    target_link_directories(nextsimlib PUBLIC ${xios_LIBRARIES})
    target_link_libraries(nextsimlib PUBLIC xios HDF5::HDF5)
    target_include_directories(
        nextsimlib
        PRIVATE ${xios_INCLUDES} ${xios_EXTERNS}/blitz/ ${xios_EXTERNS}/rapidxml/include
    )
endif()

# Do we want to include the OASIS interface? It'll only work if we have MPI
if (ENABLE_OASIS)
    if (ENABLE_MPI)
        message(STATUS "Building with OASIS support")
        find_package(OASIS REQUIRED)
        pkg_search_module(NETCDF REQUIRED netcdf-fortran)
        target_compile_definitions(nextsimlib PUBLIC USE_OASIS)

        target_link_directories(nextsimlib PUBLIC ${NETCDF_LIBRARY_DIRS})
        target_link_directories(nextsimlib PUBLIC ${OASIS_LIBRARIES})
        target_link_libraries(nextsimlib PUBLIC oasis.cbind mct mpeu psmile.MPI1 scrip netcdff netcdf)
        target_include_directories(nextsimlib PUBLIC ${OASIS_INCLUDES})
    else ()
        message(FATAL_ERROR "Cannot build with OASIS support, because MPI is not enabled" .)
    endif ()
endif ()

# OpenMP
if(WITH_THREADS)
    find_package(OpenMP REQUIRED)
    target_link_libraries(nextsimlib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Set an empty list of sources
set(NextsimSources "")
# Set an empty list of include directories
set(NextsimIncludeDirs "")
# At least one test shell script runs the full binary. This variable contains the path to it.
set(NEXTSIM_BINARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/nextsim")

# Build the core model. Defines the 'parse_modules' target
add_subdirectory("core")

# Build all components
foreach(compo ${CodeComponents})
    add_subdirectory("${compo}")
endforeach()

target_sources(nextsimlib PRIVATE ${NextsimSources})
target_include_directories(nextsimlib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${NextsimIncludeDirs}")

# main executable
add_executable(nextsim "./core/src/main.cpp")
target_link_libraries(nextsim PRIVATE nextsimlib)
# nextsimlib does not expose the needed includes so they need to be set again
target_include_directories(nextsim PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${NextsimIncludeDirs}")

# Generate the restart files that don't require additional data
add_subdirectory(run)
