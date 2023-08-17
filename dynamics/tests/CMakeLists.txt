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
