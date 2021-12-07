cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_STANDARD 11)

project(framework_dg)

find_package(netCDF REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(Catch2 REQUIRED)

# To add netCDF to a target:
# target_include_directories(target PUBLIC ${netCDF_INCLUDE_DIR})
# target_link_directories(target PUBLIC ${netCDF_LIB_DIR})
# target_link_libraries(target LINK_PUBLIC netcdf-cxx4)

# Set the location of the ipp files used by ModuleLoader for the main build
set(ModuleLoaderIppTargetDirectory
"${CMAKE_CURRENT_SOURCE_DIR}/src/modules/generated/")

# Set the list of components that define modules. Component subdirectories are
# contained in this directory level and contain src/modules/ and 
# src/modules/include/ subdirectories 
set(ModularComponents
"physics"
)

# Build the main model
add_subdirectory(src)

# Testing
if (NOT(CMAKE_BUILD_TYPE STREQUAL "Release"))
    add_subdirectory(test)
endif()
