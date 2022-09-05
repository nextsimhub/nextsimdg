# Build the main neXtSIM model

INCLUDE_DIRECTORIES(${Nextsim_INCLUDE_DIRS})
LINK_LIBRARIES(Eigen3::Eigen ${NextsimDGLib})

add_executable(example1-sasipmesh example1-sasipmesh.cpp)


