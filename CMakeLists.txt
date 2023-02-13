cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)

project(nextsim_dg)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build")

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

find_package(MPI REQUIRED)



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

#find_library(xios libxios.a REQUIRED)# HINTS ENV LD_LIBRARY_PATH LIBRARY_PATH)

include(CMakePrintHelpers)
cmake_print_variables(CMAKE_SOURCE_DIR)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#pkg_search_module(XIOS REQUIRED xios IMPORTED_TARGET)
find_package(xios REQUIRED)

cmake_print_variables(xios_INCLUDES)

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
    "${netCDF_INCLUDE_DIR}"
    "${MPI_CXX_INCLUDE_PATH}"
    "${xios_INCLUDES}"
    "${xios_EXTERNS}/blitz/"
    "${xios_EXTERNS}/rapidxml/include"
    )
target_link_directories(nextsim PUBLIC "${netCDF_LIB_DIR}" "${xios_LIBRARIES}")
target_link_libraries(nextsim LINK_PUBLIC Boost::program_options Boost::log "${NSDG_NetCDF_Library}" Eigen3::Eigen MPI::MPI_CXX xios "-lhdf5_hl -lhdf5 -lifcore")
