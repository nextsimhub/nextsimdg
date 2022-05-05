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

add_executable(testNextsimPhysics
    "NextsimPhysics_test.cpp"
    "${ModulesDir}/NextsimPhysics.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${CoreSourceDir}/ElementData.cpp"
    "${CoreSourceDir}/PrognosticElementData.cpp"
    "${ModulesDir}/HiblerConcentration.cpp"
    "${ModulesDir}/ThermoIce0ED.cpp"
    "${ModulesDir}/IPhysics1dModule.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${ModulesDir}/IIceAlbedoModule.cpp"
    "${ModulesDir}/IThermodynamicsModule.cpp"
    "${ModulesDir}/IIceOceanHeatFluxModule.cpp"
    "${ModulesDir}/IConcentrationModelModule.cpp"
    
    )
target_include_directories(testNextsimPhysics PRIVATE
    "${CoreSourceDir}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    "${netCDF_INCLUDE_DIR}"
    )
target_link_libraries(testNextsimPhysics PRIVATE "${Boost_LIBRARIES}" Catch2::Catch2)

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
    "${ModulesDir}/ThermoIce0Growth.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${ModulesDir}/IceTemperatureModule.cpp"
    "${ModulesDir}/ThermoIce0Temperature.cpp"
    "${ModulesDir}/FluxCalculationModule.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    )
target_include_directories(testIceGrowth PRIVATE
    "${CoreSourceDir}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testIceGrowth PRIVATE "${Boost_LIBRARIES}" Catch2::Catch2 Eigen3::Eigen)
    
add_executable(testTI0Temperature
    "ThermoIce0Temperature_test.cpp"
    "${CoreSourceDir}/Configurator.cpp"
    "${CoreSourceDir}/ConfiguredModule.cpp"
    "${CoreSourceDir}/ModelArray.cpp"
    "${CoreSourceDir}/ModelComponent.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${ModulesDir}/ThermoIce0Temperature.cpp"
    "${ModulesDir}/FluxCalculationModule.cpp"
    "${ModulesDir}/FiniteElementFluxes.cpp"
    )
target_include_directories(testTI0Temperature PRIVATE
    "${CoreSourceDir}"
    "${CoreModulesDir}"
    "${SourceDir}"
    "${ModulesDir}"
    )
target_link_libraries(testTI0Temperature PRIVATE "${Boost_LIBRARIES}" Catch2::Catch2 Eigen3::Eigen)

#add_executable(testFEFluxes
#    "FiniteElementFluxes_test.cpp"
#    "${CoreSourceDir}/Configurator.cpp"
#    "${CoreSourceDir}/ConfiguredModule.cpp"
#    "${CoreSourceDir}/ModelArray.cpp"
#    "${CoreSourceDir}/ModelComponent.cpp"
#    "${CoreModulesDir}/IFreezingPointModule.cpp"
#    "${ModulesDir}/FiniteElementFluxes.cpp"
#    )
#target_include_directories(testFEFluxes PRIVATE
#    "${CoreSourceDir}"
#    "${CoreModulesDir}"
#    "${SourceDir}"
#    "${ModulesDir}"
#    )
#target_link_libraries(testFEFluxes PRIVATE "${Boost_LIBRARIES}" Catch2::Catch2 Eigen3::Eigen)
