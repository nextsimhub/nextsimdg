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


add_executable(test1-sasipmesh test1-sasipmesh.cpp
 "${CoreDir}/ModelArray.cpp"
 "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
 "${SRC_DIR}/stopwatch.cpp"
 "${SRC_DIR}/Interpolations.cpp"
 "${SRC_DIR}/ParametricMesh.cpp"
 "${SRC_DIR}/DGTransport.cpp"
 "${SRC_DIR}/ParametricMap.cpp"
 "${SRC_DIR}/cgParametricMomentum.cpp"
 )
target_include_directories(test1-sasipmesh PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(test1-sasipmesh LINK_PUBLIC Eigen3::Eigen)

add_executable(test2-sasipmesh test2-sasipmesh.cpp
 "${CoreDir}/ModelArray.cpp"
 "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
 "${SRC_DIR}/stopwatch.cpp"
 "${SRC_DIR}/DGTransport.cpp"
 "${SRC_DIR}/ParametricMap.cpp"
 "${SRC_DIR}/Interpolations.cpp"
 "${SRC_DIR}/ParametricMesh.cpp"
 "${SRC_DIR}/cgParametricMomentum.cpp"
 )
target_include_directories(test2-sasipmesh PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(test2-sasipmesh LINK_PUBLIC Eigen3::Eigen)
