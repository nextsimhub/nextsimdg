cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 17)

project(nextsim_dg)

set (NETCDF_CXX "YES")
find_package(netCDF REQUIRED)
if ("${CMAKE_HOST_SYSTEM_NAME}" STREQUAL "Darwin")
    set(NSDG_NetCDF_Library "netcdf-cxx4")
else()
    set(NSDG_NetCDF_Library "netcdf_c++4")
endif()
find_package(Boost COMPONENTS program_options log REQUIRED)
find_package(Catch2 REQUIRED)
find_package(Eigen3 3.3 REQUIRED)

# To add netCDF to a target:
# target_include_directories(target PUBLIC ${netCDF_INCLUDE_DIR})
# target_link_directories(target PUBLIC ${netCDF_LIB_DIR})
# target_link_libraries(target LINK_PUBLIC "${NSDG_NetCDF_Library}")

# Set the location of the ipp files used by ModuleLoader for the main build
set(ModuleLoaderIppTargetDirectory
"${CMAKE_CURRENT_SOURCE_DIR}/core/src/modules/generated/")

file(MAKE_DIRECTORY "${ModuleLoaderIppTargetDirectory}")
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
    "${Boost_INCLUDE_DIRS}"
    "${ModuleLoaderIppTargetDirectory}"
    "${NextsimIncludeDirs}"
    "${netCDF_INCLUDE_DIR}"
    )
target_link_directories(nextsim PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(nextsim LINK_PUBLIC ${Boost_LIBRARIES} "${NSDG_NetCDF_Library}" Eigen3::Eigen)
