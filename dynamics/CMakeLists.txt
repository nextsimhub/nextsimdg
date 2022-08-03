### Main CMAKE-FILE for DG Nextsim


################################################## CMAKE options
cmake_minimum_required(VERSION 3.1.0)
project(NextsimDG CXX)


################################################## C++ Options
set(CMAKE_CXX_STANDARD 11)

#
# Setup compiler flags and compile definitions:
#
SET(NEXTSIM_CXX_FLAGS
  "-std=c++11 -pedantic -Wall -Wno-unknown-pragmas -march=native" CACHE STRING "Configuration flags"
  )

###
OPTION(WITH_THREADS      "Build with support for openmp" OFF)

find_package(OpenMP)
if (OPENMP_FOUND)
   IF(WITH_THREADS)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
   endif()
endif()



SEPARATE_ARGUMENTS(NEXTSIM_CXX_FLAGS)
ADD_COMPILE_OPTIONS(${NEXTSIM_CXX_FLAGS})

################################################## This project
### directory to sources
SET (NextsimDG_DIR ${CMAKE_CURRENT_SOURCE_DIR})

### direcetory for binaries -> per default to build
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY .)


### Name of the NextsimDG Library
set(NextsimDGLib NextsimDG)

### By standard, build a shared library
OPTION(BUILD_SHARED_LIBS       "Build shared libraries" ON)

### By standard, release
set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

if (DEFINED CMAKE_CXX_FLAGS_RELEASE_INIT
   AND "${CMAKE_CXX_FLAGS_RELEASE_INIT}" STREQUAL " ${CMAKE_CXX_FLAGS_RELEASE}")
     # Overwrite the init values choosen by CMake
#     if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
         set(CMAKE_CXX_FLAGS_RELEASE "-O3" CACHE STRING "" FORCE)
#     endif()
 endif()


################################################## Requirements of NextsimDG

### Use Eigen3
find_package (Eigen3 3.3 REQUIRED) #

### Use clang-format
OPTION(USE_CLANG_FORMAT "Formats the codebase" ON)
if(USE_CLANG_FORMAT)
  FILE(GLOB_RECURSE ALL_CXX_SOURCE_FILES src/*.cpp src/*.hpp  tests/*.cpp tests/*.hpp applications/*.cpp  applications/*.hpp)

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

################################################## Includes
LIST(APPEND Nextsim_INCLUDE_DIRS "${NextsimDG_DIR}/src/include")

################################################## Compile the Library
add_subdirectory(src)


################################################## Compile the test cases
add_subdirectory(tests)

################################################## Compile larger application examples
add_subdirectory(applications)
