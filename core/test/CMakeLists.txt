# Build the unit, integration and model tests for neXtSIM

set(INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    "../src/include"
    "../../core/src"
    "../../core/src/modules"
    "../../core/src/discontinuousgalerkin"
    # "../../core/src/include"
    )

add_executable(testIterator
    "Iterator_test.cpp"
    )
target_include_directories(testIterator PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testIterator PRIVATE nextsimlib doctest::doctest)

add_executable(testCommandLineParser
    "CommandLineParser_test.cpp"
    "ArgV.cpp"
    )
target_include_directories(testCommandLineParser PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testCommandLineParser PRIVATE nextsimlib doctest::doctest)

add_executable(testConfigurator
    "Configurator_test.cpp"
    "ArgV.cpp"
    )
target_include_directories(testConfigurator PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testConfigurator PRIVATE nextsimlib doctest::doctest)

add_executable(testConfiguredModule
    "ConfiguredModule_test.cpp"
    "ArgV.cpp"
)
target_include_directories(testConfiguredModule PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testConfiguredModule PRIVATE nextsimlib doctest::doctest)

add_executable(testTimer
    "Timer_test.cpp"
    )
target_include_directories(testTimer PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testTimer PRIVATE nextsimlib doctest::doctest)

add_executable(testScopedTimer
    "ScopedTimer_test.cpp"
    "../src/ScopedTimer.cpp"
    )
target_include_directories(testScopedTimer PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testScopedTimer PRIVATE nextsimlib doctest::doctest)

set(INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    "../src/include"
    "../../core/src"
    "../../core/src/modules"
    "./testmodelarraydetails/"
    # "../../core/src/include"
    )

add_executable(testModelArray
    "ModelArray_test.cpp"
    "testmodelarraydetails/ModelArrayDetails.cpp"
    )

target_include_directories(testModelArray PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testModelArray PRIVATE nextsimlib doctest::doctest)

add_executable(testDevGrid
    "DevGrid_test.cpp"
    )

set(INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    "../src/include"
    "../../core/src"
    "../../core/src/modules"
    "../../core/src/discontinuousgalerkin"
    "${netCDF_INCLUDE_DIR}"
    # "../../core/src/include"
    )

target_include_directories(testDevGrid PRIVATE ${INCLUDE_DIRS})
target_link_directories(testDevGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testDevGrid PRIVATE nextsimlib doctest::doctest "${NSDG_NetCDF_Library}" Eigen3::Eigen)

set(INCLUDE_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    "../src/include"
    "../../core/src"
    "../../core/src/modules"
    "../../core/src/finitevolume"
    "${netCDF_INCLUDE_DIR}"
    # "../../core/src/include"
    )

add_executable(testRectGrid
    "RectGrid_test.cpp"
    )

target_include_directories(testRectGrid PRIVATE ${INCLUDE_DIRS})
target_link_directories(testRectGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testRectGrid PRIVATE Boost::program_options doctest::doctest "${NSDG_NetCDF_Library}" Eigen3::Eigen)

add_executable(testParaGrid
    "ParaGrid_test.cpp"
    )

target_include_directories(testParaGrid PRIVATE ${INCLUDE_DIRS})
target_link_directories(testParaGrid PUBLIC "${netCDF_LIB_DIR}")
target_link_libraries(testParaGrid PRIVATE Boost::program_options doctest::doctest "${NSDG_NetCDF_Library}" Eigen3::Eigen)

add_executable(testModelComponent
    "ModelComponent_test.cpp"
)

target_include_directories(testModelComponent PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testModelComponent PRIVATE nextsimlib doctest::doctest)

add_executable(testTime
    "Time_test.cpp"
)
target_include_directories(testTime PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testTime PRIVATE nextsimlib doctest::doctest)

add_executable(testNewModelArrayRef
    "NewModelArrayRef_test.cpp"
)

target_include_directories(testNewModelArrayRef PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testNewModelArrayRef PRIVATE nextsimlib doctest::doctest)

add_executable(testPrognosticData
    "PrognosticData_test.cpp"
    )
target_include_directories(testPrognosticData PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testPrognosticData PRIVATE nextsimlib doctest::doctest)
