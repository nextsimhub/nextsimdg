# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")

set(CoreDir "../../core/src/")

include_directories(${INCLUDE_DIR})

add_executable(nextsimMEVP NextsimMevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/parametricTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/parametricTools.cpp" 
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_link_libraries(nextsimMEVP LINK_PUBLIC Eigen3::Eigen)
target_include_directories(nextsimMEVP PRIVATE "${CoreDir}" "${CoreDir}/${ModelArrayStructure}")

add_executable(benchmark_mehlmann_sasip_mevp benchmark_mehlmann_sasip_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/parametricTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/parametricTools.cpp" 
    "${SRC_DIR}/cgParametricMomentum.cpp" 
    )
target_include_directories(benchmark_mehlmann_sasip_mevp PRIVATE "${CoreDir}")

add_executable(benchmark_mehlmann_sasip_meb benchmark_mehlmann_sasip_meb.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/parametricTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/parametricTools.cpp" 
    "${SRC_DIR}/cgParametricMomentum.cpp" 
)
target_include_directories(benchmark_mehlmann_sasip_meb PRIVATE "${CoreDir}")

add_executable(benchmark_box_sasip_mevp benchmark_box_sasip_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/parametricTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/parametricTools.cpp" 
    "${SRC_DIR}/cgParametricMomentum.cpp" 
)
target_include_directories(benchmark_box_sasip_mevp PRIVATE "${CoreDir}")

target_link_libraries(benchmark_mehlmann_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_mehlmann_sasip_meb LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_box_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
