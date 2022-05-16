# Build the main neXtSIM model

INCLUDE_DIRECTORIES(${Nextsim_INCLUDE_DIRS})
LINK_LIBRARIES(Eigen3::Eigen ${NextsimDGLib})


add_executable(benchmark_box_mevp benchmark_box_mevp.cpp)
add_executable(benchmark_mehlmann_mevp benchmark_mehlmann_mevp.cpp)
add_executable(benchmark_meb benchmark_MEB_CGvel.cpp)
add_executable(benchmark_mmeb benchmark_mMEB.cpp)
add_executable(compression_meb compression_MEB_CGvel.cpp)
add_executable(compression_vero compression_Vero_MEB_CGvel.cpp)

