# Build the main neXtSIM model

INCLUDE_DIRECTORIES(${Nextsim_INCLUDE_DIRS})
LINK_LIBRARIES(Eigen3::Eigen ${NextsimDGLib})


add_executable(benchmark     benchmark.cpp)
add_executable(benchmark_box_mevp benchmark_box_mevp.cpp)
add_executable(benchmark_mehlmann_mevp benchmark_mehlmann_mevp.cpp)
add_executable(benchmark_meb benchmark_MEB_CGvel.cpp)

