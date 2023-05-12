# Build the unit, integration and model tests for neXtSIM

set(CoreSrc "../src")
set(SRC_DIR "${CoreSrc}")
set(INCLUDE_DIR "${CoreSrc}/include")
set(CoreModulesDir "${CoreSrc}/modules")
set(ModelArrayDetails "testmodelarraydetails")
set(FVModelArrayDetails "${CoreSrc}/finitevolume")

# add_executable(testexe
#   test/TestSrc.cpp
#   otherSource.cpp)
#target_link_libraries(testexe PRIVATE doctest::doctest)

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

add_executable(testTimeClasses
    "Time_test.cpp"
    "${CoreSrc}/Time.cpp"
)
target_include_directories(testTimeClasses PRIVATE "${CoreSrc}")
target_link_libraries(testTimeClasses PRIVATE doctest::doctest)

add_executable(testNewModelArrayRef
    "NewModelArrayRef_test.cpp"
    "${CoreSrc}/ModelArray.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
)

target_include_directories(testNewModelArrayRef PRIVATE "${CoreSrc}" "${CoreSrc}/${ModelArrayStructure}")
target_link_libraries(testNewModelArrayRef PRIVATE doctest::doctest Eigen3::Eigen)

add_executable(testModelArrayRefDebugging
    "ModelArrayRefDebug_test.cpp"
    "${CoreSrc}/ModelArray.cpp"
    "${FVModelArrayDetails}/ModelArrayDetails.cpp"
)
target_compile_definitions(testModelArrayRefDebugging PRIVATE DEBUG_MODELARRAYREF)
target_include_directories(testModelArrayRefDebugging PRIVATE "${CoreSrc}" "${CoreSrc}/${ModelArrayStructure}")
target_link_libraries(testModelArrayRefDebugging PRIVATE doctest::doctest Eigen3::Eigen)

# PrognosticData (and hopefully that alone) requires code from the physics tree
set(PhysSourceDir "${CoreSrc}/../../physics/src/")
set(PhysModulesDir "${CoreSrc}/../../physics/src/modules")

add_executable(testPrognosticData
    "PrognosticData_test.cpp"
    "${CoreSrc}/PrognosticData.cpp"
    "${CoreSrc}/CommonRestartMetadata.cpp"
    "${CoreSrc}/Configurator.cpp"
    "${CoreSrc}/ConfiguredModule.cpp"
    "${CoreSrc}/ModelArray.cpp"
    "${CoreSrc}/ModelComponent.cpp"
    "${CoreSrc}/ModelMetadata.cpp"
    "${CoreSrc}/MissingData.cpp"
    "${CoreSrc}/ParaGridIO.cpp"
    "${CoreModulesDir}/IFreezingPointModule.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${CoreSrc}/NZLevels.cpp"
    "${PhysSourceDir}/IceGrowth.cpp"
    "${PhysSourceDir}/IceMinima.cpp"
    "${PhysSourceDir}/SlabOcean.cpp"
    "${PhysModulesDir}/FiniteElementFluxes.cpp"
    "${PhysModulesDir}/FiniteElementSpecHum.cpp"
    "${PhysModulesDir}/FluxCalculationModule.cpp"
    "${PhysModulesDir}/OceanBoundaryModule.cpp"
    "${PhysModulesDir}/ConstantOceanBoundary.cpp"
    "${PhysModulesDir}/ConfiguredOcean.cpp"
    "${PhysModulesDir}/FluxConfiguredOcean.cpp"
    "${PhysModulesDir}/TOPAZOcean.cpp"
    "${PhysModulesDir}/AtmosphereBoundaryModule.cpp"
    "${PhysModulesDir}/ConstantAtmosphereBoundary.cpp"
    "${PhysModulesDir}/ConfiguredAtmosphere.cpp"
    "${PhysModulesDir}/FluxConfiguredAtmosphere.cpp"
    "${PhysModulesDir}/ERA5Atmosphere.cpp"
    "${PhysModulesDir}/LateralIceSpreadModule.cpp"
    "${PhysModulesDir}/HiblerSpread.cpp"
    "${PhysModulesDir}/IceThermodynamicsModule.cpp"
    "${PhysModulesDir}/ThermoIce0.cpp"
    "${PhysModulesDir}/ThermoWinton.cpp"
    "${PhysModulesDir}/IceOceanHeatFluxModule.cpp"
    "${PhysModulesDir}/BasicIceOceanHeatFlux.cpp"
    "${PhysModulesDir}/IIceAlbedoModule.cpp"
    "${PhysModulesDir}/SMUIceAlbedo.cpp"
    "${PhysModulesDir}/CCSMIceAlbedo.cpp"
    "${PhysModulesDir}/SMU2IceAlbedo.cpp"
    "${CoreSrc}/Time.cpp"
    )
target_include_directories(testPrognosticData PRIVATE
    "${CoreSrc}"
    "${CoreSrc}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${PhysSourceDir}"
    "${PhysModulesDir}"
    "${netCDF_INCLUDE_DIR}"
    )
target_link_directories(testPrognosticData PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testPrognosticData PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen "${NSDG_NetCDF_Library}")

add_executable(testConfigOutput
    "ConfigOutput_test.cpp"
    "${CoreSrc}/CommonRestartMetadata.cpp"
    "${CoreSrc}/Configurator.cpp"
    "${CoreSrc}/ConfiguredModule.cpp"
    "${CoreSrc}/Logged.cpp"
    "${CoreSrc}/ModelArray.cpp"
    "${CoreSrc}/ModelComponent.cpp"
    "${CoreSrc}/ModelMetadata.cpp"
    "${CoreSrc}/MissingData.cpp"
    "${CoreSrc}/NZLevels.cpp"
    "${CoreSrc}/DevGridIO.cpp"
    "${CoreSrc}/ParaGridIO.cpp"
    "${CoreSrc}/RectGridIO.cpp"
    "${CoreSrc}/StructureFactory.cpp"
    "${CoreSrc}/Time.cpp"
    "${CoreSrc}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    "${CoreModulesDir}/IStructureModule.cpp"
    "${CoreModulesDir}/DevGrid.cpp"
    "${CoreModulesDir}/ParametricGrid.cpp"
    "${CoreModulesDir}/RectangularGrid.cpp"
    "${CoreModulesDir}/DiagnosticOutputModule.cpp"
    "${CoreModulesDir}/SimpleOutput.cpp"
    "${CoreModulesDir}/ConfigOutput.cpp"
    )
target_include_directories(testConfigOutput PRIVATE
    "${CoreSrc}"
    "${CoreSrc}/${ModelArrayStructure}"
    "${CoreModulesDir}"
    "${netCDF_INCLUDE_DIR}"
    )
target_link_directories(testConfigOutput PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testConfigOutput PRIVATE Boost::program_options Boost::log doctest::doctest Eigen3::Eigen "${NSDG_NetCDF_Library}")
