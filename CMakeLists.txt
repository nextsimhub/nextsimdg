cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)

project(nextsim_dg)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")

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

pkg_search_module(NETCDF_FORTRAN netcdf-fortran)
if (NETCDF_FORTRAN_FOUND)
    set(NetCDF_F_Library "${NETCDF_FORTRAN_LIBRARIES}")
    set(netCDF_F_LIB_DIR "${NETCDF_FORTRAN_LIBRARY_DIRS}")
    message(STATUS "Fortran lib: ${NetCDF_F_Library} ${NETCDF_FORTRAN_LIBRARY_DIRS}")
endif()



set(OASIS_Libraries "-loasis.cbind -lpsmile.MPI1 -lscrip -lmct -lmpeu -lmpi_mpifh -lgfortran -lm")
set(OASIS_include_dir "/cluster/home/maisonnave/SASIP/oasis3-mct_5/oa3-mct_gfortran/include")
set(OASIS_lib_dir "/cluster/home/maisonnave/SASIP/oasis3-mct_5/oa3-mct_gfortran/lib")

OPTION(WITH_THREADS      "Build with support for openmp" OFF)

find_package(OpenMP)
if (OPENMP_FOUND)
   IF(WITH_THREADS)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
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

#find_package(Catch2 REQUIRED)
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/include)
add_library(Catch2::Catch2 IMPORTED INTERFACE)
set_property(TARGET Catch2::Catch2 PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${CATCH_INCLUDE_DIR}")

find_package(Eigen3 3.4 REQUIRED)

# To add netCDF to a target:
# target_include_directories(target PUBLIC ${netCDF_INCLUDE_DIR})
# target_link_directories(target PUBLIC ${netCDF_LIB_DIR})
# target_link_libraries(target LINK_PUBLIC "${NSDG_NetCDF_Library}")

# Set the list of components that define modules. Component subdirectories are
# contained in this directory level and contain src/modules/ and 
# src/modules/include/ subdirectories 
set(ModularComponents
"physics"
)

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

target_compile_options(nextsim PUBLIC "$<$<CONFIG:DEBUG>:-DDEBUG>")
target_include_directories(nextsim PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NextsimIncludeDirs}"
    "${netCDF_INCLUDE_DIR}"
    "${OASIS_include_dir}"
    )
target_link_directories(nextsim PUBLIC "${OASIS_lib_dir}" "${netCDF_LIB_DIR}" "${netCDF_F_LIB_DIR}")
target_link_libraries(nextsim LINK_PUBLIC MPI::MPI_CXX Boost::program_options Boost::log "${OASIS_Libraries}" "${NSDG_NetCDF_Library}" "${NetCDF_F_Library}" Eigen3::Eigen)
