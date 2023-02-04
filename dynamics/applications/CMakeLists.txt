# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")

include_directories(${INCLUDE_DIR})

# add_executable(benchmark_mehlmann_orcamesh_mevp benchmark_mehlmann_orcamesh_mevp.cpp
#     "${SRC_DIR}/stopwatch.cpp"
#     "${SRC_DIR}/parametricTransport.cpp"
#     "${SRC_DIR}/Interpolations.cpp"
#     "${SRC_DIR}/ParametricMesh.cpp"
#     "${SRC_DIR}/cgParametricMomentum.cpp" 
#     )
# 
add_executable(benchmark_mehlmann_sasip_mevp benchmark_mehlmann_sasip_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"
    )
add_executable(benchmark_cartesian_globe benchmark_cartesian_globe.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"
    )
add_executable(benchmark_spherical_globe benchmark_spherical_globe.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"    
    )
add_executable(test_spherical_globe test_spherical_globe.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"    
    )


  # add_executable(benchmark_mehlmann_sasip_meb benchmark_mehlmann_sasip_meb.cpp
#     "${SRC_DIR}/stopwatch.cpp"
#     "${SRC_DIR}/parametricTransport.cpp"
#     "${SRC_DIR}/Interpolations.cpp"
#     "${SRC_DIR}/ParametricMesh.cpp"
#     "${SRC_DIR}/cgParametricMomentum.cpp" 
# )
# add_executable(benchmark_box_sasip_mevp benchmark_box_sasip_mevp.cpp
#     "${SRC_DIR}/stopwatch.cpp"
#     "${SRC_DIR}/parametricTransport.cpp"
#     "${SRC_DIR}/Interpolations.cpp"
#     "${SRC_DIR}/ParametricMesh.cpp"
#     "${SRC_DIR}/cgParametricMomentum.cpp" 
# )
# 
target_link_libraries(benchmark_mehlmann_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_cartesian_globe LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_spherical_globe LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(test_spherical_globe LINK_PUBLIC Eigen3::Eigen)
#target_link_libraries(benchmark_mehlmann_orcamesh_mevp LINK_PUBLIC Eigen3::Eigen)
#target_link_libraries(benchmark_mehlmann_sasip_meb LINK_PUBLIC Eigen3::Eigen)
#target_link_libraries(benchmark_box_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
