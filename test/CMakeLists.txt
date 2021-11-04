# Build the unit, integration and model tests for neXtSIM

set(SRC_DIR "../src")
set(INCLUDE_DIR "${SRC_DIR}/include")

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

add_executable(testNextsimPhysics
    "NextsimPhysics_test.cpp"
    "${SRC_DIR}/NextsimPhysics.cpp"
    "${SRC_DIR}/PhysicsData.cpp"
    "${SRC_DIR}/ModuleLoader.cpp"
    "${SRC_DIR}/Configured.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/SMUIceAlbedo.cpp"
    "${SRC_DIR}/CCSMIceAlbedo.cpp"
    "${SRC_DIR}/SMU2IceAlbedo.cpp"
    "${SRC_DIR}/BasicIceOceanHeatFlux.cpp"
    "${SRC_DIR}/PrognosticData.cpp"
    )
target_include_directories(testNextsimPhysics PRIVATE ${ModuleLoaderIppTargetDirectory} "${SRC_DIR}")
target_link_libraries(testNextsimPhysics PRIVATE ${Boost_LIBRARIES} Catch2::Catch2)

#add_executable(testThermoIce0
#    "ThermoIce0_test.cpp"
#    "${SRC_DIR}/ThermoIce0.cpp"
#    "${SRC_DIR}/PrognosticData.cpp"
#    )
#target_link_libraries(testThermoIce0 PRIVATE Catch2::Catch2)

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


