# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")

set(CoreDir "../../core/src/")

include_directories(${INCLUDE_DIR})

add_executable(createNextsimMeshData createNextsimMeshData.cpp
    "${CoreDir}/modules/ParametricGrid.cpp"
    "${CoreDir}/modules/RectangularGrid.cpp"
    "${CoreDir}/modules/DevGrid.cpp"
    "${CoreDir}/modules/IStructureModule.cpp"
    "${CoreDir}/CommonRestartMetadata.cpp"
    "${CoreDir}/Configurator.cpp"
    "${CoreDir}/ParaGridIO.cpp"
    "${CoreDir}/MissingData.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/ModelMetadata.cpp"
    "${CoreDir}/NZLevels.cpp"
    "${CoreDir}/Time.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_link_libraries(createNextsimMeshData LINK_PUBLIC Eigen3::Eigen Boost::program_options "${NSDG_NetCDF_Library}")
target_link_directories(createNextsimMeshData PUBLIC "${netCDF_LIB_DIR}")
target_include_directories(createNextsimMeshData PRIVATE "${CoreDir}" "${CoreDir}/modules" "${CoreDir}/${ModelArrayStructure}" "${netCDF_INCLUDE_DIR}")

add_executable(nextsimMEVP NextsimMevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${CoreDir}/modules/ParametricGrid.cpp"
    "${CoreDir}/modules/RectangularGrid.cpp"
    "${CoreDir}/modules/DevGrid.cpp"
    "${CoreDir}/modules/IStructureModule.cpp"
    "${CoreDir}/CommonRestartMetadata.cpp"
    "${CoreDir}/Configurator.cpp"
    "${CoreDir}/MissingData.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/ModelMetadata.cpp"
    "${CoreDir}/NZLevels.cpp"
    "${CoreDir}/ParaGridIO.cpp"
    "${CoreDir}/Time.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_link_libraries(nextsimMEVP LINK_PUBLIC Eigen3::Eigen Boost::program_options "${NSDG_NetCDF_Library}")
target_link_directories(nextsimMEVP PUBLIC "${netCDF_LIB_DIR}")
target_include_directories(nextsimMEVP PRIVATE "${CoreDir}" "${CoreDir}/modules" "${CoreDir}/${ModelArrayStructure}" "${netCDF_INCLUDE_DIR}")

add_executable(benchmark_freedrift_meb benchmark_freedrift_meb.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    )
#add_executable(benchmark_freedrift_mevp benchmark_freedrift_mevp.cpp
#    "${SRC_DIR}/stopwatch.cpp"
#    "${SRC_DIR}/parametricTransport.cpp"
#    "${SRC_DIR}/Interpolations.cpp"
#    "${SRC_DIR}/ParametricMesh.cpp"
#    "${SRC_DIR}/ParametricMap.cpp"
#    "${SRC_DIR}/cgParametricMomentum.cpp"
#    )

add_executable(benchmark_mehlmann_sasip_mevp benchmark_mehlmann_sasip_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"
    )
add_executable(benchmark_mehlmann_orcamesh_mevp benchmark_mehlmann_orcamesh_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"
    )
add_executable(benchmark_mehlmann_orcamesh_mevp_spherical benchmark_mehlmann_orcamesh_mevp_spherical.cpp
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


add_executable(benchmark_mehlmann_sasip_meb benchmark_mehlmann_sasip_meb.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"
    )
add_executable(benchmark_box_sasip_mevp benchmark_box_sasip_mevp.cpp
    "${SRC_DIR}/stopwatch.cpp"
    "${SRC_DIR}/DGTransport.cpp"
    "${SRC_DIR}/Interpolations.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${SRC_DIR}/ParametricMap.cpp"
    "${SRC_DIR}/cgParametricMomentum.cpp"
    "${SRC_DIR}/VectorManipulations.cpp"
    )

target_link_libraries(benchmark_mehlmann_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_cartesian_globe LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_spherical_globe LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(test_spherical_globe LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_mehlmann_orcamesh_mevp LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_mehlmann_orcamesh_mevp_spherical LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_mehlmann_sasip_meb LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(benchmark_box_sasip_mevp LINK_PUBLIC Eigen3::Eigen)
