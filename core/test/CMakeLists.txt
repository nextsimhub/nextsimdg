# Build the unit, integration and model tests for neXtSIM

set(CoreSrc "../src")
set(SRC_DIR "${CoreSrc}")
set(INCLUDE_DIR "${CoreSrc}/include")
set(CoreModulesDir "${CoreSrc}/modules")
set(ModelArrayDetails "testmodelarraydetails")

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
    "${SRC_DIR}/Time.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/Logged.cpp"
    "${SRC_DIR}/Configurator.cpp"
    
    )
target_include_directories(testIterator PUBLIC "${SRC_DIR}")
target_link_libraries(testIterator PRIVATE doctest::doctest Boost::program_options Boost::log)

add_executable(testCommandLineParser
    "CommandLineParser_test.cpp"
    "ArgV.cpp"
    "${SRC_DIR}/CommandLineParser.cpp"
    "${SRC_DIR}/ConfigurationHelpPrinter.cpp"
    )
target_include_directories(testCommandLineParser PUBLIC "${SRC_DIR}")
target_link_libraries(testCommandLineParser LINK_PUBLIC Boost::program_options Boost::log doctest::doctest)

add_executable(testConfigurator
    "Configurator_test.cpp"
    "ArgV.cpp"
    "${SRC_DIR}/Configurator.cpp"
    )
target_include_directories(testConfigurator PUBLIC "${SRC_DIR}")
target_link_libraries(testConfigurator LINK_PUBLIC Boost::program_options Boost::log doctest::doctest)

add_executable(testConfiguredModule
    "ConfiguredModule_test.cpp"
    "ArgV.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ConfiguredModule.cpp"
)
target_include_directories(testConfiguredModule PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${SRC_DIR}")
target_link_libraries(testConfiguredModule PRIVATE doctest::doctest Boost::program_options Boost::log)

add_executable(testTimer
    "Timer_test.cpp"
    "${SRC_DIR}/Timer.cpp"
    )
target_link_libraries(testTimer PRIVATE doctest::doctest)
target_include_directories(testTimer PRIVATE "${SRC_DIR}")

add_executable(testScopedTimer
    "ScopedTimer_test.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/ScopedTimer.cpp"
    )
target_link_libraries(testScopedTimer PRIVATE doctest::doctest)
target_include_directories(testScopedTimer PRIVATE "${SRC_DIR}")

add_executable(testModelArray
    "ModelArray_test.cpp"
    "${CoreSrc}/ModelArray.cpp"
    "${ModelArrayDetails}/ModelArrayDetails.cpp"
    )

target_include_directories(testModelArray PRIVATE "${CoreSrc}" "${ModelArrayDetails}")
target_link_libraries(testModelArray PRIVATE doctest::doctest Eigen3::Eigen)

add_executable(testDevGrid
    "DevGrid_test.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${SRC_DIR}/CommonRestartMetadata.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ConfiguredModule.cpp"
    "${SRC_DIR}/DevGridIO.cpp"
    "${SRC_DIR}/MissingData.cpp"
    "${SRC_DIR}/ModelArray.cpp"
    "${SRC_DIR}/ModelMetadata.cpp"
    "${SRC_DIR}/NZLevels.cpp"
    "${SRC_DIR}/Time.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )

target_include_directories(testDevGrid PUBLIC "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}" "${netCDF_INCLUDE_DIR}" "${CoreSrc}/${ModelArrayStructure}")
target_link_directories(testDevGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testDevGrid LINK_PUBLIC Boost::program_options Boost::log doctest::doctest "${NSDG_NetCDF_Library}" Eigen3::Eigen)

add_executable(testRectGrid
    "RectGrid_test.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${SRC_DIR}/CommonRestartMetadata.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ConfiguredModule.cpp"
    "${SRC_DIR}/RectGridIO.cpp"
    "${SRC_DIR}/MissingData.cpp"
    "${SRC_DIR}/ModelArray.cpp"
    "${SRC_DIR}/ModelMetadata.cpp"
    "${SRC_DIR}/NZLevels.cpp"
    "${SRC_DIR}/Time.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )

target_include_directories(testRectGrid PUBLIC "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}" "${netCDF_INCLUDE_DIR}" "${CoreSrc}/${ModelArrayStructure}")
target_link_directories(testRectGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testRectGrid LINK_PUBLIC Boost::program_options doctest::doctest "${NSDG_NetCDF_Library}" Eigen3::Eigen)

add_executable(testParaGrid
    "ParaGrid_test.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${SRC_DIR}/CommonRestartMetadata.cpp"
    "${SRC_DIR}/Configurator.cpp"
    "${SRC_DIR}/ConfiguredModule.cpp"
    "${SRC_DIR}/ParaGridIO.cpp"
    "${SRC_DIR}/MissingData.cpp"
    "${SRC_DIR}/ModelArray.cpp"
    "${SRC_DIR}/ModelMetadata.cpp"
    "${SRC_DIR}/NZLevels.cpp"
    "${SRC_DIR}/Time.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )

target_include_directories(testParaGrid PUBLIC "${SRC_DIR}" "${CoreModulesDir}" "${PhysicsDir}" "${PhysicsModulesDir}" "${netCDF_INCLUDE_DIR}" "${CoreSrc}/${ModelArrayStructure}")
target_link_directories(testParaGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testParaGrid LINK_PUBLIC Boost::program_options doctest::doctest "${NSDG_NetCDF_Library}" Eigen3::Eigen)

add_executable(testModelComponent
    "ModelComponent_test.cpp"
    "${CoreSrc}/ModelComponent.cpp"
    "${CoreSrc}/MissingData.cpp"
    "${CoreSrc}/Configurator.cpp"
    "${CoreSrc}/ModelArray.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
)

target_include_directories(testModelComponent PRIVATE "${CoreSrc}" "${CoreSrc}/${ModelArrayStructure}")
target_link_libraries(testModelComponent PRIVATE Boost::program_options doctest::doctest Eigen3::Eigen)

add_executable(testTime
    "Time_test.cpp"
    "${CoreSrc}/Time.cpp"
)
target_include_directories(testTime PRIVATE "${CoreSrc}")
target_link_libraries(testTime PRIVATE doctest::doctest)

add_executable(testNewModelArrayRef
    "NewModelArrayRef_test.cpp"
    "${CoreSrc}/ModelArray.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
)

target_include_directories(testNewModelArrayRef PRIVATE "${CoreSrc}" "${CoreSrc}/${ModelArrayStructure}")
target_link_libraries(testNewModelArrayRef PRIVATE doctest::doctest Eigen3::Eigen)
