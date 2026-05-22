cmake_minimum_required(VERSION 3.12)

# set globally because we have a number of independent targets
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(nextsim_dg)

# user definable build options
option(WITH_THREADS "Build with support for openmp" OFF)
option(WITH_KOKKOS "Build the Kokkos accelerated dynamical core" OFF)
option(ENABLE_MPI "Enable distributed-memory parallelization with MPI" OFF)
option(ENABLE_XIOS "Enable XIOS library for IO" OFF)
option(BUILD_TESTS "Build the tests" ON)
set(DynamicsType
    "DG2"
    CACHE STRING "Set the discretization order for the dynamics, either DG1 or DG2")
option(USE_SINGLE_PRECISION "Use float instead of double everywhere" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# all nextsim code is collected in a library to allow for reuse in the tests the code is added after
# processing the subdirectories
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
    find_package(
        NetCDF
        COMPONENTS C CXX
        REQUIRED)
    set(NSDG_NetCDF_Library "${NetCDF_LIBRARIES}")
    set(netCDF_INCLUDE_DIR "${NetCDF_INCLUDE_DIRS}")
endif()

target_include_directories(nextsimlib PUBLIC "${netCDF_INCLUDE_DIR}")
target_link_directories(nextsimlib PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(nextsimlib PUBLIC "${NSDG_NetCDF_Library}")

if(BUILD_TESTS)
    # Add the doctest header library
    set(DOCTEST_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
    add_library(doctest::doctest IMPORTED INTERFACE)
    set_property(TARGET doctest::doctest PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                  "${DOCTEST_INCLUDE_DIR}")
    enable_testing() # Enable ctest
endif()

# Eigen
find_package(Eigen3 3.4) # Try to find Eigen 3.4
if(NOT TARGET Eigen3::Eigen)
    # If Eigen 3.4 was not found, require any Eigen version after 3.4.1 and before 6
    find_package(Eigen3 3.4...5 REQUIRED NO_MODULE)
endif()
target_link_libraries(nextsimlib PUBLIC Eigen3::Eigen)

if(DEFINED PYTHON_EXECUTABLE)
    set(Python_EXECUTABLE
        PYTHON_EXECUTABLE
        CACHE STRING "")
else()
    set(Python_FIND_VIRTUALENV FIRST)
    find_package(Python COMPONENTS Interpreter)
endif()

# The location of the module_builder.py scripts
set(ScriptDirectory "${PROJECT_SOURCE_DIR}/scripts")
set(ModuleBuilderScript "${ScriptDirectory}/module_builder.py")
set(ModuleHeaderScript "${ScriptDirectory}/module_header.py")

# Dynamics type. Defaults to DG1. Available options are: DG1 Discontinuous Galerkin, degree 1 (3
# components) DG2 Discontinuous Galerkin, degree 2 (6 components)
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
target_compile_definitions(nextsimlib PRIVATE DGCOMP=${DGComp} DGSTRESSCOMP=${DGStressComp}
                                              CGDEGREE=${CGDegree})

if(USE_SINGLE_PRECISION)
    if(ENABLE_MPI)
        message(SEND_ERROR "Trying to use both MPI and single precision together, which is currently not supported. \
                Disable either ENABLE_MPI or USE_SINGLE_PRECISION.")
    endif()
    if(ENABLE_XIOS)
        message(SEND_ERROR "Trying to use both XIOS and single precision together, which is currently not supported. \
                Disable either ENABLE_XIOS or USE_SINGLE_PRECISION.")
    endif()
    if (NOT WITH_KOKKOS)
        message(WARNING "Single precision has only been tested with the Kokkos implementations.")
    endif()

    # EIGEN_FAST_MATH enables an alternative sqrt implementation when used with float.
    # This has not been sufficiently tested and may currently cause issues in the dynamics.
    target_compile_definitions(nextsimlib PRIVATE EIGEN_FAST_MATH=0)

    set(FloatType "float")
else()
    set(FloatType "double")
endif()
# set the floating point type
target_compile_definitions(nextsimlib PUBLIC "NEXTSIM_FLOAT_TYPE=${FloatType}")
message("FP type is set to ${FloatType}")

# boost
find_package(
    Boost
    COMPONENTS program_options log
    REQUIRED)
target_link_libraries(nextsimlib PUBLIC Boost::program_options Boost::log)
# Regarding Boost.Log, if our application consists of more than one modules that use it, we must
# link to the shared version. If we have a single executable or a single module that works, we may
# use the static version. By default, it is assumed that the library is built in static mode. Use
# the following definition to indicate that the code will be linked against dynamically loaded boost
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
        PUBLIC ${xios_INCLUDES} ${xios_EXTERNS}/blitz ${xios_EXTERNS}/boost_extraction/include
                ${xios_EXTERNS}/rapidxml/include
    )
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(FORTRAN_RUNTIME_LIB "-lgfortran -lmpi_mpifh")
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(FORTRAN_RUNTIME_LIB "-lifcore")
    endif()
    target_link_libraries(nextsimlib PUBLIC "${FORTRAN_RUNTIME_LIB}")
endif()

# OpenMP
if(WITH_THREADS)
    find_package(OpenMP REQUIRED)
    target_link_libraries(nextsimlib PUBLIC OpenMP::OpenMP_CXX)
endif()

# Kokkos
if(WITH_KOKKOS)
    if(ENABLE_MPI)
        message(SEND_ERROR "Trying to use both MPI and Kokkos which is currently not supported. \
                Disable either ENABLE_MPI or WITH_KOKKOS.")
    endif()
    # If nextsimlib is a shared library than kokkos need to be build with position independent code.
    get_target_property(target_type nextsimlib TYPE)
    if(${target_type} STREQUAL "SHARED_LIBRARY")
        # target specific property is difficult here set_target_properties(kokkos PROPERTIES
        # POSITION_INDEPENDENT_CODE ON)
        set(CMAKE_POSITION_INDEPENDENT_CODE True)
    endif()

    find_package(Kokkos 4.3 QUIET CONFIG)
    if(Kokkos_FOUND)
        message(STATUS "Found Kokkos: ${Kokkos_DIR} (version \"${Kokkos_VERSION}\")")
    else()
        message("No kokkos installation found. Looking for local version in \"lib/\".")
        add_subdirectory("lib/kokkos")
    endif()

    target_link_libraries(nextsimlib PUBLIC Kokkos::kokkos)
    # define to use for conditional code inclusion
    target_compile_definitions(nextsimlib PUBLIC USE_KOKKOS)

    if(Kokkos_ENABLE_CUDA OR Kokkos_ENABLE_ROCM)
        # Eigen's intrinsic based implementations are not valid on GPU
        target_compile_definitions(nextsimlib PRIVATE "EIGEN_DONT_VECTORIZE")
        # 64 bit ints are only emulated so it is better to use a smaller type
        target_compile_definitions(nextsimlib PUBLIC EIGEN_DEFAULT_DENSE_INDEX_TYPE=int)
    endif()
    if(Kokkos_ENABLE_CUDA)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
            message("Clang is currently not supported with CUDA.")
        else()
            # Disable warning caused by Eigen: "calling a __host__ function from a __host__
            # __device__ function is not allowed" It is not sufficient to disable this for just CUDA
            # code, i.e. with the generator expression: $<$<COMPILE_LANGUAGE:CUDA>:-Xcudafe
            # --diag_suppress=20014> "#68-D: integer conversion resulted in a change of sign" This
            # warning is also caused by CPU code (e.g. Slice.hpp) and should be fixed at some point.
            # "#611-D: overloaded virtual function "[...]" is only partially overridden in class
            # "[...]" Caused by child classes of ModelComponent.
            target_compile_options(
                nextsimlib
                PUBLIC "SHELL:-Xcudafe --diag_suppress=20014" "SHELL:-Xcudafe --diag_suppress=68"
                       "SHELL:-Xcudafe --diag_suppress=611")
        endif()
    endif()
    # stdandard is not set correctly for amdclang based on just the target property does not matter
    # for versions >=17 because c++17 is the default
    if(Kokkos_ENABLE_ROCM AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        target_compile_options(nextsimlib PUBLIC -std=c++17)
    endif()
    # reset to use defaults again
    unset(CMAKE_POSITION_INDEPENDENT_CODE)
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
