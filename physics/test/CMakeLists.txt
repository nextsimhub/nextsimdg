# Build the unit, integration and model tests for neXtSIM

set(CoreSourceDir "${PROJECT_SOURCE_DIR}/core/src")
set(CoreModulesDir "${CoreSourceDir}/modules")
set(SourceDir "../src")
set(ModulesDir "${SourceDir}/modules")

# add_executable(testexe
#   test/TestSrc.cpp
#   otherSource.cpp)
#target_link_libraries(testexe PRIVATE Catch2::Catch2)

# include(CTest)
# include(Catch)
# catch_discover_tests(testexe)

add_executable(testIceGrowth
    "IceGrowth_test.cpp"
    "${SourceDir}/IceGrowth.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${CoreSourceDir}/NZLevels.cpp"
    "${ModulesDir}/IceThermodynamicsModule.cpp"
    "${ModulesDir}/LateralIceSpreadModule.cpp"
    "${ModulesDir}/HiblerSpread.cpp"
    "${ModulesDir}/ThermoIce0.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testIceGrowth PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testIceGrowth PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen)

add_executable(testFEFluxes
    "FiniteElementFluxes_test.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testFEFluxes PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testFEFluxes PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen)


add_executable(testERA5Atm
    "ERA5Atm_test.cpp"
    "${ModulesDir}/ERA5Atmosphere.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/CommonRestartMetadata.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/ModelMetadata.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreSourceDir}/NZLevels.cpp"
    "${CoreSourceDir}/ParaGridIO.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    "${ModulesDir}/FluxCalculationModule.cpp"
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_compile_definitions(testERA5Atm PRIVATE TEST_FILE_SOURCE=${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(testERA5Atm PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_directories(testERA5Atm PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testERA5Atm PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen "${NSDG_NetCDF_Library}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/era5_test128x128.nc"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(testTOPAZOcn
    "TOPAZOcn_test.cpp"
    "${ModulesDir}/TOPAZOcean.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/CommonRestartMetadata.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/ModelMetadata.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreSourceDir}/NZLevels.cpp"
    "${CoreSourceDir}/ParaGridIO.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${CoreSourceDir}/Time.cpp"
    "${ModulesDir}/IceOceanHeatFluxModule.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    )
target_compile_definitions(testTOPAZOcn PRIVATE TEST_FILE_SOURCE=${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(testTOPAZOcn PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_directories(testTOPAZOcn PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testTOPAZOcn PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen "${NSDG_NetCDF_Library}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/topaz_test128x128.nc"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(testBIOHFluxes
    "BasicIceOceanFlux_test.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/Time.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_include_directories(testBIOHFluxes PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testBIOHFluxes PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen)

add_executable(testSpecHum
    "SpecificHumidity_test.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    )
target_include_directories(testSpecHum PRIVATE
    "${ModulesDir}"
    )
target_link_libraries(testSpecHum PRIVATE Catch2::Catch2)
