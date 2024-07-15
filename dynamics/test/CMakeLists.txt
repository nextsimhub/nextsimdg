# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")
set(CoreDir "../../core/src/")

include_directories(${INCLUDE_DIR})

add_executable(testDGModelArray
    "DGModelArray_test.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_include_directories(testDGModelArray PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(testDGModelArray LINK_PUBLIC doctest::doctest Eigen3::Eigen)

add_executable(testCGModelArray
    "CGModelArray_test.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_include_directories(testCGModelArray PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(testCGModelArray LINK_PUBLIC doctest::doctest Eigen3::Eigen)

add_executable(testParametricMesh
    "ParametricMesh_test.cpp"
    "FakeSmeshData.cpp"
    "${SRC_DIR}/ParametricMesh.cpp"
    "${CoreDir}/ModelArray.cpp"
    "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
target_include_directories(testParametricMesh PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(testParametricMesh LINK_PUBLIC doctest::doctest)


add_executable(testAdvection Advection_test.cpp
 "${SRC_DIR}/Interpolations.cpp"
 "${SRC_DIR}/DGTransport.cpp"
 "${SRC_DIR}/ParametricMap.cpp"
 )
target_include_directories(testAdvection PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(testAdvection LINK_PUBLIC doctest::doctest Eigen3::Eigen)

add_executable(testAdvectionPeriodicBC AdvectionPeriodicBC_test.cpp
 "${SRC_DIR}/DGTransport.cpp"
 "${SRC_DIR}/ParametricMap.cpp"
 "${SRC_DIR}/Interpolations.cpp"
 )
target_include_directories(testAdvectionPeriodicBC PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}")
target_link_libraries(testAdvectionPeriodicBC LINK_PUBLIC doctest::doctest Eigen3::Eigen)
