# Build the unit, integration and model tests for neXtSIM

set(CoreSourceDir "${PROJECT_SOURCE_DIR}/core/src")
set(CoreModulesDir "${CoreSourceDir}/modules")
set(SourceDir "../src")
set(ModulesDir "${SourceDir}/modules")

# add_executable(testexe
#   test/TestSrc.cpp
#   otherSource.cpp)
#target_link_libraries(testexe PRIVATE doctest::doctest)

add_executable(testIceGrowth
    "IceGrowth_test.cpp"
    "${SourceDir}/IceGrowth.cpp"
    "${SourceDir}/IceMinima.cpp"
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
    "${ModulesDir}/ThermoWinton.cpp"
    "${ModulesDir}/UniformOcean.cpp"
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
target_link_libraries(testIceGrowth PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)

add_executable(testThermoWinton
    "ThermoWintonTemperature_test.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreSourceDir}/NZLevels.cpp"
    "${SourceDir}/IceMinima.cpp"
    "${ModulesDir}/LateralIceSpreadModule.cpp"
    "${ModulesDir}/HiblerSpread.cpp"
    "${ModulesDir}/ThermoWinton.cpp"
    "${ModulesDir}/UniformOcean.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testThermoWinton PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testThermoWinton PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)

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
    "${ModulesDir}/MU71Albedo.cpp"
    "${ModulesDir}/UniformOcean.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testFEFluxes PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testFEFluxes PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)


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
    "${CoreSourceDir}/FileCallbackCloser.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    "${ModulesDir}/FluxCalculationModule.cpp"
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${ModulesDir}/MU71Albedo.cpp"
    "${ModulesDir}/UniformOcean.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_compile_definitions(testERA5Atm PRIVATE TEST_FILE_SOURCE=${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(testERA5Atm PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    "${netCDF_INCLUDE_DIR}"
    )
target_link_directories(testERA5Atm PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testERA5Atm PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen "${NSDG_NetCDF_Library}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/era5_test128x128.nc"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(testTOPAZOcn
    "TOPAZOcn_test.cpp"
    "${ModulesDir}/TOPAZOcean.cpp"
    "${SourceDir}/SlabOcean.cpp"
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
    "${CoreSourceDir}/FileCallbackCloser.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
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
    "${netCDF_INCLUDE_DIR}"
    )
target_link_directories(testTOPAZOcn PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testTOPAZOcn PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen "${NSDG_NetCDF_Library}")
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/topaz_test128x128.nc"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(testBIOHFluxes
    "BasicIceOceanFlux_test.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${ModulesDir}/UniformOcean.cpp"
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
target_link_libraries(testBIOHFluxes PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)

add_executable(testSpecHum
    "SpecificHumidity_test.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    )
target_include_directories(testSpecHum PRIVATE
    "${ModulesDir}"
    )
target_link_libraries(testSpecHum PRIVATE doctest::doctest)

add_executable(testUniformOcean
    "UniformOcean_test.cpp"
    "${ModulesDir}/UniformOcean.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
)
target_include_directories(testUniformOcean PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${ModulesDir}"
    )
target_link_libraries(testUniformOcean PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)

add_executable(testConstantOcn
    "ConstantOceanBoundary_test.cpp"
    "${ModulesDir}/ConstantOceanBoundary.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    "${ModulesDir}/FluxCalculationModule.cpp"
    "${ModulesDir}/IceOceanHeatFluxModule.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${ModulesDir}/MU71Albedo.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testConstantOcn PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testConstantOcn PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)

add_executable(testSlabOcn
    "SlabOcean_test.cpp"
    "${SourceDir}/SlabOcean.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testSlabOcn PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testSlabOcn PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)

add_executable(testThermoIce0
    "ThermoIce0_test.cpp"
    "${ModulesDir}/ThermoIce0.cpp"
    "${SourceDir}/IceMinima.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreSourceDir}/NZLevels.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreSourceDir}/Time.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_include_directories(testThermoIce0 PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testThermoIce0 PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)

add_executable(testIceMinima
    "IceMinima_test.cpp"
    "${SourceDir}/IceMinima.cpp"
    )
target_include_directories(testIceMinima PRIVATE "${SourceDir}")
target_link_libraries(testIceMinima PRIVATE doctest::doctest)

add_executable(testBenchmarkBoundaries
    "BenchmarkBoundaries_test.cpp"
    "${ModulesDir}/BenchmarkAtmosphere.cpp"
    "${ModulesDir}/BenchmarkOcean.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreSourceDir}/MissingData.cpp"
    "${CoreSourceDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${SourceDir}/BenchmarkCoordinates.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testBenchmarkBoundaries PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testBenchmarkBoundaries PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen)
