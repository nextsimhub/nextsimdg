# Build the unit, integration and model tests for neXtSIM

set(SRC_DIR "../src")
set(INCLUDE_DIR "${SRC_DIR}/include")
set(CoreModulesDir "${SRC_DIR}/modules")

# add_executable(testexe
#   test/TestSrc.cpp
#   otherSource.cpp)
#target_link_libraries(testexe PRIVATE Catch2::Catch2)

# include(CTest)
# include(Catch)
# catch_discover_tests(testexe)

include_directories(${INCLUDE_DIR})

add_executable(testIterator
    "Iterator_test.cpp"
    "${SRC_DIR}/Iterator.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/Logged.cpp"
    )
target_link_libraries(testIterator PRIVATE Catch2::Catch2)

add_executable(testSimpleIter
    "SimpleIterant_test.cpp"
    "${SRC_DIR}/SimpleIterant.cpp"
    "${SRC_DIR}/Iterator.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/Logged.cpp"
    )
target_link_libraries(testSimpleIter PRIVATE Catch2::Catch2)

add_executable(testCommandLineParser
    "CommandLineParser_test.cpp"
    "ArgV.cpp"
    "${SRC_DIR}/CommandLineParser.cpp"
    )
target_include_directories(testCommandLineParser PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(testCommandLineParser LINK_PUBLIC ${Boost_LIBRARIES} Catch2::Catch2)

add_executable(testConfigurator
    "Configurator_test.cpp"
    "ArgV.cpp"
    "${SRC_DIR}/Configurator.cpp"
    )
target_include_directories(testConfigurator PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(testConfigurator LINK_PUBLIC ${Boost_LIBRARIES} Catch2::Catch2)

add_executable(testEnumWrapper
    "EnumWrapper_test.cpp"
    )
target_include_directories(testEnumWrapper PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(testEnumWrapper LINK_PUBLIC ${Boost_LIBRARIES} Catch2::Catch2)

# Set the location of the test module loader classes
set(TEST_IPP_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ModuleLoaderTestModules")

add_executable(testModuleLoader
    "ModuleLoader_test.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
)
target_include_directories(testModuleLoader PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${INCLUDE_DIR}" "${TEST_IPP_DIR}")
target_link_libraries(testModuleLoader LINK_PUBLIC Catch2::Catch2)

add_executable(testConfiguredModule
    "ConfiguredModule_test.cpp"
    "ArgV.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ConfiguredModule.cpp"
)
target_include_directories(testConfiguredModule PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${INCLUDE_DIR}" "${TEST_IPP_DIR}" ${Boost_INCLUDE_DIRS})
target_link_libraries(testConfiguredModule PRIVATE Catch2::Catch2 ${Boost_LIBRARIES})

add_executable(testTimer
    "Timer_test.cpp"
    "${SRC_DIR}/Timer.cpp"
    )
target_link_libraries(testTimer PRIVATE Catch2::Catch2)
target_include_directories(testTimer PRIVATE "${SRC_DIR}")

add_executable(testScopedTimer
    "ScopedTimer_test.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/ScopedTimer.cpp"
    )
target_link_libraries(testScopedTimer PRIVATE Catch2::Catch2)
target_include_directories(testScopedTimer PRIVATE "${SRC_DIR}")

add_executable(testPrognosticData
    "PrognosticData_test.cpp"
    "${SRC_DIR}/PrognosticData.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    )
    # Set the location of the test module loader classes
set(PDTestIppDir "${CMAKE_CURRENT_SOURCE_DIR}/PrognosticDataTestModules")
target_link_libraries(testPrognosticData PRIVATE Catch2::Catch2)
target_include_directories(testPrognosticData PRIVATE "${SRC_DIR}" "${CoreModulesDir}" "${PDTestIppDir}")

set(PhysicsDir "${PROJECT_SOURCE_DIR}/physics/src")
set(PhysicsModulesDir "${PhysicsDir}/modules")
add_executable(testElementData
    "ElementData_test.cpp"
    "${SRC_DIR}/ElementData.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ConfiguredModule.cpp"
    "${SRC_DIR}/PrognosticData.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${PhysicsModulesDir}/NextsimPhysics.cpp"
    "${PhysicsModulesDir}/SMUIceAlbedo.cpp"
    "${PhysicsModulesDir}/CCSMIceAlbedo.cpp"
    "${PhysicsModulesDir}/SMU2IceAlbedo.cpp"
    "${PhysicsModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${PhysicsModulesDir}/HiblerConcentration.cpp"
    "${PhysicsModulesDir}/ThermoIce0.cpp"
    )

target_include_directories(testElementData PRIVATE "${ModuleLoaderIppTargetDirectory}" "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}")
target_link_libraries(testElementData PRIVATE "${Boost_LIBRARIES}" Catch2::Catch2)

add_executable(exampleDevGridOutput
    "DevGrid_example.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    "${SRC_DIR}/ElementData.cpp"
    "${SRC_DIR}/PrognosticData.cpp"
    "${SRC_DIR}/DevGridIO.cpp"
    "${PhysicsModulesDir}/NextsimPhysics.cpp"
    "${PhysicsModulesDir}/SMUIceAlbedo.cpp"
    "${PhysicsModulesDir}/CCSMIceAlbedo.cpp"
    "${PhysicsModulesDir}/SMU2IceAlbedo.cpp"
    "${PhysicsModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${PhysicsModulesDir}/HiblerConcentration.cpp"
    "${PhysicsModulesDir}/ThermoIce0.cpp"
    )

target_include_directories(exampleDevGridOutput PUBLIC "${ModuleLoaderIppTargetDirectory}" "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}" "${netCDF_INCLUDE_DIR}")
target_link_directories(exampleDevGridOutput PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(exampleDevGridOutput LINK_PUBLIC "${Boost_LIBRARIES}" Catch2::Catch2 "${NSDG_NetCDF_Library}")

add_executable(testDevGrid
    "DevGrid_test.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    "${SRC_DIR}/ElementData.cpp"
    "${SRC_DIR}/PrognosticData.cpp"
    "${SRC_DIR}/DevGridIO.cpp"
    "${PhysicsModulesDir}/NextsimPhysics.cpp"
    "${PhysicsModulesDir}/SMUIceAlbedo.cpp"
    "${PhysicsModulesDir}/CCSMIceAlbedo.cpp"
    "${PhysicsModulesDir}/SMU2IceAlbedo.cpp"
    "${PhysicsModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${PhysicsModulesDir}/HiblerConcentration.cpp"
    "${PhysicsModulesDir}/ThermoIce0.cpp"
    )

target_include_directories(testDevGrid PUBLIC "${ModuleLoaderIppTargetDirectory}" "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}" "${netCDF_INCLUDE_DIR}")
target_link_directories(testDevGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testDevGrid LINK_PUBLIC "${Boost_LIBRARIES}" Catch2::Catch2 "${NSDG_NetCDF_Library}")

add_executable(testRectGrid
    "RectGrid_test.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    "${SRC_DIR}/ElementData.cpp"
    "${SRC_DIR}/PrognosticData.cpp"
    "${SRC_DIR}/RectGridIO.cpp"
    "${PhysicsModulesDir}/NextsimPhysics.cpp"
    "${PhysicsModulesDir}/SMUIceAlbedo.cpp"
    "${PhysicsModulesDir}/CCSMIceAlbedo.cpp"
    "${PhysicsModulesDir}/SMU2IceAlbedo.cpp"
    "${PhysicsModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${PhysicsModulesDir}/HiblerConcentration.cpp"
    "${PhysicsModulesDir}/ThermoIce0.cpp"
    )

target_include_directories(testRectGrid PUBLIC "${ModuleLoaderIppTargetDirectory}" "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}" "${netCDF_INCLUDE_DIR}")
target_link_directories(testRectGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testRectGrid LINK_PUBLIC Boost::program_options Catch2::Catch2 "${NSDG_NetCDF_Library}")

add_executable(testStructureFactory
    "StructureFactory_test.cpp"
    "${SRC_DIR}/StructureFactory.cpp"
    "${SRC_DIR}/ElementData.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ConfiguredModule.cpp"
    "${SRC_DIR}/PrognosticData.cpp"
    "${SRC_DIR}/DevGridIO.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${PhysicsModulesDir}/NextsimPhysics.cpp"
    "${PhysicsModulesDir}/SMUIceAlbedo.cpp"
    "${PhysicsModulesDir}/CCSMIceAlbedo.cpp"
    "${PhysicsModulesDir}/SMU2IceAlbedo.cpp"
    "${PhysicsModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${PhysicsModulesDir}/HiblerConcentration.cpp"
    "${PhysicsModulesDir}/ThermoIce0.cpp"
    )

target_include_directories(testStructureFactory PRIVATE "${ModuleLoaderIppTargetDirectory}" "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}" "${netCDF_INCLUDE_DIR}")
target_link_directories(testStructureFactory PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testStructureFactory PRIVATE "${Boost_LIBRARIES}" Catch2::Catch2 "${NSDG_NetCDF_Library}")
