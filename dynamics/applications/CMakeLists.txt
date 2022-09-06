# Build the main neXtSIM model

INCLUDE_DIRECTORIES(${Nextsim_INCLUDE_DIRS})
LINK_LIBRARIES(Eigen3::Eigen ${NextsimDGLib})

add_executable(benchmark_mehlmann_sasip_mevp benchmark_mehlmann_sasip_mevp.cpp)
add_executable(benchmark_mehlmann_sasip_meb benchmark_mehlmann_sasip_meb.cpp)
add_executable(benchmark_box_sasip_mevp benchmark_box_sasip_mevp.cpp)



