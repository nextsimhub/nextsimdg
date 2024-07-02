cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)

project(nextsim_dg)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(PkgConfig)
pkg_search_module(NETCDF_CXX4 netcdf-cxx4)
if (NETCDF_CXX4_FOUND)
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

OPTION(WITH_THREADS      "Build with support for openmp" OFF)

find_package(OpenMP)
if (OPENMP_FOUND)
    IF(WITH_THREADS)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

option(ENABLE_XIOS "Enable XIOS library for IO" OFF)
if (ENABLE_XIOS)
    message(STATUS "Building with XIOS library")
    find_package(xios REQUIRED QUIET)
    if (XIOS_FOUND)
        find_package(HDF5 REQUIRED COMPONENTS C HL)
    endif()
endif()

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

# Add the doctest header library
set(DOCTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
add_library(doctest::doctest IMPORTED INTERFACE)
set_property(TARGET doctest::doctest PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${DOCTEST_INCLUDE_DIR}")

find_package(Eigen3 3.4 REQUIRED)

option(ENABLE_MPI "Enable distributed-memory parallelization with MPI" OFF)
if(ENABLE_MPI)
    find_package(MPI REQUIRED COMPONENTS C CXX)
endif()

# if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
#     set(FORTRAN_RUNTIME_LIB "-lgfortran")
# elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
#     set(FORTRAN_RUNTIME_LIB "-lifcore")
# endif()

set(Python_FIND_VIRTUALENV FIRST)
find_package (Python COMPONENTS Interpreter)

# To add netCDF to a target:
# target_include_directories(target PUBLIC ${netCDF_INCLUDE_DIR})
# target_link_directories(target PUBLIC ${netCDF_LIB_DIR})
# target_link_libraries(target LINK_PUBLIC "${NSDG_NetCDF_Library}")

# Set the list of components that define modules. Component subdirectories are
# contained in this directory level and contain src/modules/ and
# src/modules/include/ subdirectories
#set(ModularComponents
#"physics"
#)

# The location of the module_builder.py scripts
set(ScriptDirectory "${PROJECT_SOURCE_DIR}/scripts")
set(ModuleBuilderScript "${ScriptDirectory}/module_builder.py")
set(ModuleHeaderScript "${ScriptDirectory}/module_header.py")

# Dynamics type. Defaults to DG1.
# Available options are:
# DG1 Discontinuous Galerkin, degree 1 (3 components)
# DG2 Discontinuous Galerkin, degree 2 (6 components)
set(isDG FALSE)

# Replace a missing dynamics type with DG1
if (NOT DynamicsType)
    message("Defaulting dynamics to DG1")
    set(DynamicsType "DG1")
endif()

message("Dynamics is set to ${DynamicsType}")

if (DynamicsType STREQUAL "DG1")
    set(isDG TRUE)
    set(DGComp 3)
    set(CGDegree 2)
endif()

if (DynamicsType STREQUAL "DG2")
    set(isDG TRUE)
    set(DGComp 6)
    set(CGDegree 2)
endif()

# Set the components which provide source or object code to the main model
set(CodeComponents
    "physics"
    )

if (isDG)
# Add the DG dynamics subdirectory to the list of code components
    set(CodeComponents
        "${CodeComponents}"
        "dynamics"
        )

# Set the number of DG stress components given the CG degree
    if (CGDegree EQUAL 1)
        set(DGStressComp 3)
    elseif(CGDegree EQUAL 2)
        set(DGStressComp 8)
    else()
        message("Invalid value of CGDegree. Valid values are 1–2")
    endif()

# The exact selection of dimensions and Types available to ModelArray
    set(ModelArrayStructure "discontinuousgalerkin")
endif()

# Set an empty list of sources
set(NextsimSources "")

# Set an empty list of include directories
set(NextsimIncludeDirs "")

# At least one test shell script runs the full binary. This variable contains the path to it.
set(NEXTSIM_BINARY_PATH "${CMAKE_CURRENT_BINARY_DIR}/nextsim")

# Build the core model. Defines the 'parse_modules' target
add_subdirectory(core)


# Build all components
foreach(compo ${CodeComponents})
    add_subdirectory("${compo}")
endforeach()

add_library(nextsimlib SHARED ${NextsimSources})
target_compile_definitions(nextsimlib
    PRIVATE
    $<$<BOOL:${ENABLE_MPI}>:USE_MPI>
    $<$<BOOL:${ENABLE_XIOS}>:USE_XIOS>
    DGCOMP=${DGComp}
    DGSTRESSCOMP=${DGStressComp}
    CGDEGREE=${CGDegree}
    )
target_include_directories(nextsimlib
    PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NextsimIncludeDirs}"
    $<$<BOOL:${ENABLE_XIOS}>:${xios_INCLUDES}>
    $<$<BOOL:${ENABLE_XIOS}>:${xios_EXTERNS}/blitz/>
    $<$<BOOL:${ENABLE_XIOS}>:${xios_EXTERNS}/rapidxml/include>
    PUBLIC
    "${netCDF_INCLUDE_DIR}"
    )

target_link_directories(nextsimlib PUBLIC "${netCDF_LIB_DIR}" $<$<BOOL:${ENABLE_MPI}>:${xios_LIBRARIES}>)
target_link_libraries(nextsimlib PUBLIC
    $<$<BOOL:${ENABLE_XIOS}>:xios>
    $<$<BOOL:${ENABLE_XIOS}>:HDF5::HDF5>
    Boost::program_options Boost::log "${NSDG_NetCDF_Library}" Eigen3::Eigen
    $<$<BOOL:${ENABLE_MPI}>:MPI::MPI_C> $<$<BOOL:${ENABLE_MPI}>:MPI::MPI_CXX>
    "${FORTRAN_RUNTIME_LIB}"
    )

add_executable(nextsim ./core/src/main.cpp)
target_compile_definitions(nextsim
    PRIVATE
    $<$<BOOL:${ENABLE_MPI}>:USE_MPI>
    $<$<BOOL:${ENABLE_XIOS}>:USE_XIOS>
    )
target_include_directories(nextsim PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NextsimIncludeDirs}"
    )
target_link_libraries(nextsim PUBLIC nextsimlib)

# Generate the restart files that don't require additional data
add_subdirectory(run)
