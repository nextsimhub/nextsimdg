cmake_minimum_required(VERSION 3.12)

project(nextsim_dg)

OPTION(WITH_THREADS "Build with support for openmp" OFF)
OPTION(WITH_KOKKOS "Build the Kokkos accelerated dynamical core" OFF)
option(ENABLE_MPI "Enable distributed-memory parallelization with MPI" OFF)
option(ENABLE_XIOS "Enable XIOS library for IO" OFF)

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

# Set the components which provide source or object code to the main model
set(CodeComponents
    "physics"
    "dynamics"
    )

# The exact selection of dimensions and Types available to ModelArray
if(True)
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
target_compile_features(nextsimlib PUBLIC cxx_std_17)
target_compile_definitions(nextsimlib
    PRIVATE
    $<$<BOOL:${ENABLE_MPI}>:USE_MPI>
    $<$<BOOL:${ENABLE_XIOS}>:USE_XIOS>
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

if (WITH_THREADS)
    find_package(OpenMP REQUIRED)
    target_link_libraries(nextsimlib PUBLIC OpenMP::OpenMP_CXX)
endif()

IF (WITH_KOKKOS)
    #if nextsimlib is a shared library than kokkos need to be build
    # with position independent code
    get_target_property(target_type nextsimlib TYPE)
    if (${target_type} STREQUAL "SHARED_LIBRARY")
        # target specific property is difficult here
        #set_target_properties(kokkos PROPERTIES POSITION_INDEPENDENT_CODE ON)
        set(CMAKE_POSITION_INDEPENDENT_CODE True)
    endif()

    find_package(Kokkos QUIET)
    IF (NOT Kokkos_FOUND)
        message("No kokkos installation found. Looking for local version in \"lib/\".")
        add_subdirectory("lib/kokkos")
    endif ()
    target_link_libraries(nextsimlib PUBLIC Kokkos::kokkos)
    
    if (Kokkos_ENABLE_CUDA)
        target_compile_definitions(nextsimlib PRIVATE "EIGEN_DONT_VECTORIZE")
        # disable warning caused by eigen
        # "calling a __host__ function from a __host__ __device__ function is not allowed"
        target_compile_options(nextsimlib PUBLIC -Xcudafe --diag_suppress=20014)#-w
    endif()
    # reset to use defaults again
    unset(CMAKE_POSITION_INDEPENDENT_CODE)
endif()

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
