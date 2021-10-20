cmake_minimum_required(VERSION 3.1.0)

set(CMAKE_CXX_STANDARD 11)

project(framework_dg)

find_package(netCDF REQUIRED)
find_package(Boost REQUIRED)
find_package(Catch2 REQUIRED)

# To add netCDF to a target:
# target_include_directories(target PUBLIC ${netCDF_INCLUDE_DIR})
# target_link_directories(target PUBLIC ${netCDF_LIB_DIR})
# target_link_libraries(target LINK_PUBLIC netcdf-cxx4)

# Build the main model
add_subdirectory(src)

# Testing
if (NOT(CMAKE_BUILD_TYPE STREQUAL "Release"))
	add_subdirectory(test)
endif()
