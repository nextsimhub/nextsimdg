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

include_directories("../src/include")

add_executable(example1-sasipmesh "example1-sasipmesh.cpp")
target_link_libraries(example1-sasipmesh PRIVATE dynamicslib)

add_executable(example2-sasipmesh "example2-sasipmesh.cpp")
target_link_libraries(example2-sasipmesh PRIVATE dynamicslib)

add_executable(example3-orcamesh "example3-orcamesh.cpp")
target_link_libraries(example3-orcamesh  PRIVATE dynamicslib)

add_executable(example4-orcamesh-spherical "example4-orcamesh-spherical.cpp")
target_link_libraries(example4-orcamesh-spherical  PRIVATE dynamicslib)

add_executable(example4-orcamesh-spherical-greenland "example4-orcamesh-spherical-greenland.cpp")
target_link_libraries(example4-orcamesh-spherical-greenland  PRIVATE dynamicslib)
