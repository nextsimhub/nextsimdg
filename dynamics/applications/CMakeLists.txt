set(COMMON_INCLUDE_DIRS
    "../../core/src/discontinuousgalerkin"
    "../../core/src"
    "../../core/src/modules"
    "../src"
    "../src/include"
    )

add_executable(createNextsimMeshData "createNextsimMeshData.cpp")
target_include_directories(createNextsimMeshData PRIVATE ${COMMON_INCLUDE_DIRS} ${netCDF_INCLUDE_DIR})
target_link_libraries(createNextsimMeshData LINK_PUBLIC nextsimlib)

set(sources
    "../src/stopwatch.cpp"
    "../src/DGTransport.cpp"
    "../src/Interpolations.cpp"
    "../src/ParametricMesh.cpp"
    "../src/ParametricMap.cpp"
    "../src/cgParametricMomentum.cpp"
    "../src/VectorManipulations.cpp"
    )

include_directories("../src/include")

add_library(dynamicslib SHARED ${sources})

add_executable(benchmark_mehlmann_sasip_mevp "benchmark_mehlmann_sasip_mevp.cpp")
target_link_libraries(benchmark_mehlmann_sasip_mevp PRIVATE dynamicslib)

add_executable(benchmark_mehlmann_orcamesh_mevp "benchmark_mehlmann_orcamesh_mevp.cpp")
target_link_libraries(benchmark_mehlmann_orcamesh_mevp PRIVATE dynamicslib)

add_executable(benchmark_mehlmann_orcamesh_mevp_spherical "benchmark_mehlmann_orcamesh_mevp_spherical.cpp")
target_link_libraries(benchmark_mehlmann_orcamesh_mevp_spherical PRIVATE dynamicslib)

add_executable(benchmark_cartesian_globe "benchmark_cartesian_globe.cpp")
target_link_libraries(benchmark_cartesian_globe PRIVATE dynamicslib)

add_executable(benchmark_spherical_globe "benchmark_spherical_globe.cpp")
target_link_libraries(benchmark_spherical_globe PRIVATE dynamicslib)

add_executable(test_spherical_globe "test_spherical_globe.cpp")
target_link_libraries(test_spherical_globe PRIVATE dynamicslib)

add_executable(benchmark_mehlmann_sasip_meb "benchmark_mehlmann_sasip_meb.cpp")
target_link_libraries(benchmark_mehlmann_sasip_meb PRIVATE dynamicslib)

add_executable(benchmark_box_sasip_mevp "benchmark_box_sasip_mevp.cpp")
target_link_libraries(benchmark_box_sasip_mevp PRIVATE dynamicslib)
