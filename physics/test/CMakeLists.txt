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

add_executable(testConfiguredAtm
    "ConfiguredAtmosphere_test.cpp"
    "${ModulesDir}/ConfiguredAtmosphere.cpp"
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
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testConfiguredAtm PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testConfiguredAtm PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen)

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
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testConstantOcn PRIVATE
    "${CoreSourceDir}"
    "${CoreSourceDir}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testConstantOcn PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen)

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
target_link_libraries(testSlabOcn PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen)
