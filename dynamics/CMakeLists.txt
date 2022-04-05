### Main CMAKE-FILE for DG Nextsim

cmake_minimum_required(VERSION 3.1.0)
# Project definition
set(PROJECT_NAME NextsimDG)
project(${PROJECT_NAME} CXX)

# CMake options
enable_testing()
set(CMAKE_BUILD_TYPE RelWithDebInfo)
# Set RPATH behaviour
set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
# Define RPATH for executables via a relative expression to enable a
# fully relocatable package
if(APPLE)
set(basePoint @loader_path)
else()
set(basePoint $ORIGIN)
endif()
include(GNUInstallDirs)
file(RELATIVE_PATH relDir
 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
 ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
# On most Unix-based platforms the $ORIGIN placeholder represents
# the location of the binary in which the RPATH is embedded
set(CMAKE_INSTALL_RPATH ${basePoint}/${relDir})

# C++ options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Setup compiler flags and compile definitions
add_compile_options(-Wpedantic -Wall -Wno-unknown-pragmas)

# Required packages
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Optional features
option(ENABLE_OPENMP       "Use OpenMP for shared-memory parallelization" OFF)
option(ENABLE_CLANG_FORMAT "Formats the codebase" ON)

if(ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
  if(OpenMP_CXX_FOUND)
    if("${OpenMP_CXX_LIBRARIES}" STREQUAL "")
      message(STATUS "No OpenMP library needs to be linked against")
    else()
      message(STATUS "Will link against OpenMP libraries: ${OpenMP_CXX_LIBRARIES}")
    endif()
  else()
    message(WARNING "OpenMP was requested, but not found.")
  endif()
endif()

if(ENABLE_CLANG_FORMAT)
  file(GLOB_RECURSE ALL_CXX_SOURCE_FILES src/*.cpp src/*.hpp tests/*.cpp tests/*.hpp applications/*.cpp  applications/*.hpp)

  # Adding clang-format target if executable is found
  find_program(CLANG_FORMAT "clang-format")
  if(CLANG_FORMAT)
    add_custom_target(
      clang-format ALL
      COMMAND ${CLANG_FORMAT}
      -i
      -style=WebKit
      ${ALL_CXX_SOURCE_FILES}
      )
  else()
    message(STATUS "Could not find clang-format")
  endif()
endif()

### Name of the NextsimDG library
set(NextsimDGLib NextsimDG)

# Compile the library
add_subdirectory(src)
# Compile larger application examples
add_subdirectory(applications)
# Compile the test cases
add_subdirectory(tests)
