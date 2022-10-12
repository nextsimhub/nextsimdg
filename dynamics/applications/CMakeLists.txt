# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")

include_directories(${INCLUDE_DIR})

add_executable(benchmark_mehlmann_sasip_mevp benchmark_mehlmann_sasip_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/parametricTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/parametricTools.cpp" 
    "${SRC_DIR}/cgParametricMomentum.cpp" 
    )
add_executable(benchmark_mehlmann_sasip_meb benchmark_mehlmann_sasip_meb.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/parametricTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/parametricTools.cpp" 
    "${SRC_DIR}/cgParametricMomentum.cpp" 
)
add_executable(benchmark_box_sasip_mevp benchmark_box_sasip_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/parametricTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/parametricTools.cpp" 
    "${SRC_DIR}/cgParametricMomentum.cpp" 
)

target_link_libraries(benchmark_mehlmann_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_mehlmann_sasip_meb LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_box_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
