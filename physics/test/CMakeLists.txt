# Build the unit, integration and model tests for neXtSIM

set(COMMON_INCLUDE_DIRS
    "../../core/src/discontinuousgalerkin"
    "../../core/src"
    "../../core/src/modules"
    "../src"
    "../src/modules"
    )

add_executable(testIceGrowth "IceGrowth_test.cpp")
target_include_directories(testIceGrowth PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testIceGrowth PRIVATE nextsimlib doctest::doctest)

add_executable(testThermoWinton "ThermoWintonTemperature_test.cpp")
target_include_directories(testThermoWinton PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testThermoWinton PRIVATE nextsimlib doctest::doctest)

add_executable(testFEFluxes "FiniteElementFluxes_test.cpp")
target_include_directories(testFEFluxes PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testFEFluxes PRIVATE nextsimlib doctest::doctest)

add_executable(testConfiguredAtm "ConfiguredAtmosphere_test.cpp")
target_include_directories(testConfiguredAtm PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testConfiguredAtm PRIVATE nextsimlib doctest::doctest)

add_executable(testBIOHFluxes "BasicIceOceanFlux_test.cpp")
target_include_directories(testBIOHFluxes PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testBIOHFluxes PRIVATE nextsimlib doctest::doctest)

add_executable(testSpecHum "SpecificHumidity_test.cpp")
target_include_directories(testSpecHum PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testSpecHum PRIVATE nextsimlib doctest::doctest)

add_executable(testConstantOcn "ConstantOceanBoundary_test.cpp")
target_include_directories(testConstantOcn PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testConstantOcn PRIVATE nextsimlib doctest::doctest)

add_executable(testSlabOcn "SlabOcean_test.cpp")
target_include_directories(testSlabOcn PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testSlabOcn PRIVATE nextsimlib doctest::doctest)

add_executable(testIceMinima "IceMinima_test.cpp")
target_include_directories(testIceMinima PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(testIceMinima PRIVATE nextsimlib doctest::doctest)
