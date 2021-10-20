# Build the main neXtSIM model

set (INC "${CMAKE_CURRENT_SOURCE_DIR}/include")


### Set the list of all cpp files
FILE(GLOB Nextsim_Sources *.cpp)


###
INCLUDE_DIRECTORIES(${Nextsim_INCLUDE_DIRS})

add_library(${NextsimDGLib} ${Nextsim_Sources})

target_link_libraries (${NextsimDGLib} Eigen3::Eigen)



