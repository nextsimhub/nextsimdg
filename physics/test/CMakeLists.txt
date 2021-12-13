# Build the unit, integration and model tests for neXtSIM

set(CommonSourceDir "${PROJECT_SOURCE_DIR}/core/src")
set(CommonIncludeDir "${SourceDir}/include")
set(SourceDir "../src")
set(ModulesDir "${SourceDir}/modules")

# add_executable(testexe
#   test/TestSrc.cpp
#   otherSource.cpp)
#target_link_libraries(testexe PRIVATE Catch2::Catch2)

# include(CTest)
# include(Catch)
# catch_discover_tests(testexe)

include_directories("${CommonSourceDir}" "${SourceDir}" "${ModulesDir}")

add_executable(testNextsimPhysics
    "NextsimPhysics_test.cpp"
    "${ModulesDir}/NextsimPhysics.cpp"
    "${CommonSourceDir}/ModuleLoader.cpp"
    "${CommonSourceDir}/Configurator.cpp"
    "${CommonSourceDir}/ConfiguredModule.cpp"
    "${ModulesDir}/SMUIceAlbedo.cpp"
    "${ModulesDir}/CCSMIceAlbedo.cpp"
    "${ModulesDir}/SMU2IceAlbedo.cpp"
    "${ModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${CommonSourceDir}/PrognosticData.cpp"
    "${ModulesDir}/HiblerConcentration.cpp"
    "${ModulesDir}/ThermoIce0.cpp"
    )
target_include_directories(testNextsimPhysics PRIVATE ${ModuleLoaderIppTargetDirectory} "${SourceDir}")
target_link_libraries(testNextsimPhysics PRIVATE ${Boost_LIBRARIES} Catch2::Catch2)

#add_executable(testThermoIce0
#    "ThermoIce0_test.cpp"
#    "${SourceDir}/ThermoIce0.cpp"
#    "${SourceDir}/PrognosticData.cpp"
#    )
#target_link_libraries(testThermoIce0 PRIVATE Catch2::Catch2)
