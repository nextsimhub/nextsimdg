link_libraries(${NextsimDGLib})

add_executable(benchmark_box_mevp benchmark_box_mevp.cpp)
add_executable(benchmark_mehlmann_mevp benchmark_mehlmann_mevp.cpp)
add_executable(benchmark_meb benchmark_MEB_CGvel.cpp)

install(TARGETS benchmark_box_mevp benchmark_mehlmann_mevp benchmark_meb
  DESTINATION ${CMAKE_INSTALL_BINDIR})
