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
    "${ModulesDir}/IceThermodynamicsModule.cpp"
    "${ModulesDir}/LateralIceSpreadModule.cpp"
    "${ModulesDir}/HiblerSpread.cpp"
    "${ModulesDir}/ThermoIce0.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${ModulesDir}/FluxCalculationModule.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${ModulesDir}/IceOceanHeatFluxModule.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${SourceDir}/AtmosphereOceanState.cpp"
    "${ModulesDir}/AtmosphereStateModule.cpp"
    "${ModulesDir}/AtmosphereState.cpp"
    "${ModulesDir}/DummyAtmosphereState.cpp"
    "${ModulesDir}/OceanStateModule.cpp"
    "${ModulesDir}/OceanState.cpp"
    "${ModulesDir}/DummyOceanState.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testIceGrowth PRIVATE
    "${CoreSourceDir}"
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
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${ModulesDir}/IceOceanHeatFluxModule.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${SourceDir}/AtmosphereOceanState.cpp"
    "${ModulesDir}/AtmosphereStateModule.cpp"
    "${ModulesDir}/AtmosphereState.cpp"
    "${ModulesDir}/DummyAtmosphereState.cpp"
    "${ModulesDir}/OceanStateModule.cpp"
    "${ModulesDir}/OceanState.cpp"
    "${ModulesDir}/DummyOceanState.cpp"
    "${CoreSourceDir}/Time.cpp"
    )
target_include_directories(testFEFluxes PRIVATE
    "${CoreSourceDir}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testFEFluxes PRIVATE Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen)

add_executable(testSpecHum
    "SpecificHumidity_test.cpp"
    "${ModulesDir}/FiniteElementSpecHum.cpp"
    )
target_include_directories(testSpecHum PRIVATE
    "${ModulesDir}"
    )
target_link_libraries(testSpecHum PRIVATE Catch2::Catch2)
