set(COMMON_INCLUDE_DIRS
    "../../core/src/discontinuousgalerkin"
    "../../core/src"
    "../../core/src/modules"
    "../src"
    "../src/include"
    )

add_executable(dgma_test "DGModelArray_test.cpp")
target_include_directories(dgma_test PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(dgma_test PRIVATE nextsimlib doctest::doctest)

add_executable(cgma_test "CGModelArray_test.cpp")
target_include_directories(cgma_test PRIVATE ${COMMON_INCLUDE_DIRS})
target_link_libraries(cgma_test PRIVATE nextsimlib doctest::doctest)

add_executable(testParametricMesh "ParametricMesh_test.cpp" "FakeSmeshData.cpp")
target_include_directories(testParametricMesh PRIVATE ${COMMON_INCLUDE_DIRS} ".")
target_link_libraries(testParametricMesh PRIVATE nextsimlib doctest::doctest)
target_compile_definitions(testParametricMesh PRIVATE TEST_FILE_SOURCE=${CMAKE_CURRENT_BINARY_DIR})
file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/25km_NH.smesh"
    DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
