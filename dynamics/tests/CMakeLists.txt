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


add_executable(example1-sasipmesh example1-sasipmesh.cpp
 "${CoreDir}/ModelArray.cpp"
 "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
 "${SRC_DIR}/stopwatch.cpp"
 "${SRC_DIR}/Interpolations.cpp"
 "${SRC_DIR}/ParametricMesh.cpp"
 "${SRC_DIR}/DGTransport.cpp"
 "${SRC_DIR}/ParametricMap.cpp"
 "${SRC_DIR}/cgParametricMomentum.cpp"
 )
target_include_directories(example1-sasipmesh PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(example1-sasipmesh LINK_PUBLIC Eigen3::Eigen)

add_executable(example2-sasipmesh example2-sasipmesh.cpp
 "${CoreDir}/ModelArray.cpp"
 "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
 "${SRC_DIR}/stopwatch.cpp"
 "${SRC_DIR}/DGTransport.cpp"
 "${SRC_DIR}/ParametricMap.cpp"
 "${SRC_DIR}/Interpolations.cpp"
 "${SRC_DIR}/ParametricMesh.cpp"
 "${SRC_DIR}/cgParametricMomentum.cpp"
 )
target_include_directories(example2-sasipmesh PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(example2-sasipmesh LINK_PUBLIC Eigen3::Eigen)

#add_executable(example3-orcamesh example3-orcamesh.cpp
#"${SRC_DIR}/stopwatch.cpp"
# "${SRC_DIR}/DGTransport.cpp"
# "${SRC_DIR}/ParametricMap.cpp"
#"${SRC_DIR}/Interpolations.cpp"
#"${SRC_DIR}/ParametricMesh.cpp"
#"${SRC_DIR}/cgParametricMomentum.cpp"
#)
#target_include_directories(example3-orcamesh PRIVATE "${CoreDir}")
#target_link_libraries(example3-orcamesh LINK_PUBLIC Eigen3::Eigen)
#
#add_executable(example4-orcamesh-spherical example4-orcamesh-spherical.cpp
#"${SRC_DIR}/stopwatch.cpp"
#"${SRC_DIR}/DGTransport.cpp"
#"${SRC_DIR}/Interpolations.cpp"
#"${SRC_DIR}/ParametricMesh.cpp"
#"${SRC_DIR}/ParametricMap.cpp"
#"${SRC_DIR}/cgParametricMomentum.cpp"
#)
#
#add_executable(example4-orcamesh-spherical-greenland example4-orcamesh-spherical-greenland.cpp
#"${SRC_DIR}/stopwatch.cpp"
#"${SRC_DIR}/DGTransport.cpp"
#"${SRC_DIR}/Interpolations.cpp"
#"${SRC_DIR}/ParametricMesh.cpp"
#"${SRC_DIR}/ParametricMap.cpp"
#"${SRC_DIR}/cgParametricMomentum.cpp"
#)
#
#target_link_libraries(example1-sasipmesh LINK_PUBLIC Eigen3::Eigen)
#target_link_libraries(example2-sasipmesh LINK_PUBLIC Eigen3::Eigen)
#target_link_libraries(example3-orcamesh  LINK_PUBLIC Eigen3::Eigen)
#target_link_libraries(example4-orcamesh-spherical  LINK_PUBLIC Eigen3::Eigen)
#target_link_libraries(example4-orcamesh-spherical-greenland  LINK_PUBLIC Eigen3::Eigen)
