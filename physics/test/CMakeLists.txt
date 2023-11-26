set(COMMON_INCLUDE_DIRS
    "../../core/src/discontinuousgalerkin"
    "../../core/src"
    "../../core/src/modules"
    "../src"
    "../src/modules"
    )

set(ModulesRoot "../src/modules")

include_directories(${COMMON_INCLUDE_DIRS})

add_executable(testIceGrowth "IceGrowth_test.cpp")
target_link_libraries(testIceGrowth PRIVATE nextsimlib doctest::doctest)

add_executable(testThermoWinton "ThermoWintonTemperature_test.cpp")
target_include_directories(testThermoWinton PRIVATE "${ModulesRoot}/IceThermodynamicsModule")
target_link_libraries(testThermoWinton PRIVATE nextsimlib doctest::doctest)

add_executable(testFEFluxes "FiniteElementFluxes_test.cpp")
target_include_directories(testFEFluxes PRIVATE "${ModulesRoot}/FluxCalculationModule")
target_link_libraries(testFEFluxes PRIVATE nextsimlib doctest::doctest)

add_executable(testBIOHFluxes "BasicIceOceanFlux_test.cpp")
target_include_directories(testBIOHFluxes PRIVATE "${ModulesRoot}/IceOceanHeatFluxModule")
target_link_libraries(testBIOHFluxes PRIVATE nextsimlib doctest::doctest)

add_executable(testSpecHum "SpecificHumidity_test.cpp")
target_include_directories(testSpecHum PRIVATE "${ModulesRoot}/FluxCalculationModule")
target_link_libraries(testSpecHum PRIVATE nextsimlib doctest::doctest)

add_executable(testConstantOcn "ConstantOceanBoundary_test.cpp")
target_link_libraries(testConstantOcn PRIVATE nextsimlib doctest::doctest)

add_executable(testSlabOcn "SlabOcean_test.cpp")
target_link_libraries(testSlabOcn PRIVATE nextsimlib doctest::doctest)

add_executable(testThermoIce0 "ThermoIce0_test.cpp")
target_include_directories(testThermoIce0 PRIVATE "${ModulesRoot}/IceThermodynamicsModule")
target_link_libraries(testThermoIce0 PRIVATE nextsimlib doctest::doctest)

add_executable(testIceMinima "IceMinima_test.cpp")
target_link_libraries(testIceMinima PRIVATE nextsimlib doctest::doctest)

add_executable(testERA5Atm "ERA5Atm_test.cpp")
target_include_directories(testERA5Atm PRIVATE "${ModulesRoot}/AtmosphereBoundaryModule")
target_compile_definitions(testERA5Atm PRIVATE TEST_FILE_SOURCE=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(testERA5Atm PRIVATE nextsimlib doctest::doctest)
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/era5_test128x128.nc"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(testTOPAZOcn "TOPAZOcn_test.cpp")
target_compile_definitions(testTOPAZOcn PRIVATE TEST_FILE_SOURCE=${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(testTOPAZOcn PRIVATE nextsimlib doctest::doctest)
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/topaz_test128x128.nc"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(testUniformOcean "UniformOcean_test.cpp")
target_link_libraries(testUniformOcean PRIVATE nextsimlib doctest::doctest)

add_executable(testBenchmarkBoundaries "BenchmarkBoundaries_test.cpp")
target_include_directories(testBenchmarkBoundaries PRIVATE "${ModulesRoot}/AtmosphereBoundaryModule")
target_link_libraries(testBenchmarkBoundaries PRIVATE nextsimlib doctest::doctest)
