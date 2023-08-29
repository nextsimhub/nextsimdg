# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")
set(CoreDir "../../core/src/")

include_directories(${INCLUDE_DIR})

add_executable(dgma_test
    "DGModelArray_test.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_include_directories(dgma_test PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(dgma_test LINK_PUBLIC doctest::doctest Eigen3::Eigen)

add_executable(cgma_test
    "CGModelArray_test.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_include_directories(cgma_test PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(cgma_test LINK_PUBLIC doctest::doctest Eigen3::Eigen)
