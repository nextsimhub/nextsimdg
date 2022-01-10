# Build the main neXtSIM model

INCLUDE_DIRECTORIES(${Nextsim_INCLUDE_DIRS})
LINK_LIBRARIES(Eigen3::Eigen ${NextsimDGLib})


add_executable(benchmark     benchmark.cpp)
add_executable(benchmark_evp evp.cpp)


