# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")
set(CoreDir "../../core/src/")

include_directories(${INCLUDE_DIR})

if(NOT ENABLE_MPI)
    add_executable(
        testDGModelArray
        "DGModelArray_test.cpp"
        "${CoreDir}/ModelArray.cpp"
        "${CoreDir}/ModelArraySlice.cpp"
        "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
    target_include_directories(
        testDGModelArray
        PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}"
    )
    target_link_libraries(testDGModelArray LINK_PUBLIC doctest::doctest Eigen3::Eigen)

    add_executable(
        testDGVectorHolder
        "dgVectorHolder_test.cpp"
        "${CoreDir}/ModelArray.cpp"
        "${CoreDir}/ModelArraySlice.cpp"
        "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
    target_include_directories(
        testDGVectorHolder
        PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}"
    )
    target_link_libraries(testDGVectorHolder LINK_PUBLIC doctest::doctest Eigen3::Eigen)

    add_executable(
        testCGModelArray
        "CGModelArray_test.cpp"
        "${CoreDir}/ModelArray.cpp"
        "${CoreDir}/ModelArraySlice.cpp"
        "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
    target_include_directories(
        testCGModelArray
        PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}"
    )
    target_link_libraries(testCGModelArray LINK_PUBLIC doctest::doctest Eigen3::Eigen)

    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/25km_NH.smesh" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

    add_executable(
        testParametricMesh
        "ParametricMesh_test.cpp"
        "FakeSmeshData.cpp"
        "${SRC_DIR}/ParametricMesh.cpp"
        "${CoreDir}/ModelArray.cpp"
        "${CoreDir}/ModelArraySlice.cpp"
        "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
    )
    target_include_directories(
        testParametricMesh
        PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}"
    )
    target_link_libraries(testParametricMesh LINK_PUBLIC doctest::doctest Eigen3::Eigen)

    add_executable(
        testParametricMeshArea
        "ParametricMeshArea_test.cpp"
        "${SRC_DIR}/ParametricMesh.cpp"
        "${CoreDir}/ModelArray.cpp"
        "${CoreDir}/${ModelArrayStructure}/ModelArrayDetails.cpp"
        "${CoreDir}/ModelArraySlice.cpp"
    )
    target_include_directories(
        testParametricMeshArea
        PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}"
    )
    target_link_libraries(testParametricMeshArea LINK_PUBLIC doctest::doctest Eigen3::Eigen)

    add_executable(
        testAdvection
        Advection_test.cpp
        "${SRC_DIR}/Interpolations.cpp"
        "${SRC_DIR}/DGTransport.cpp"
        "${SRC_DIR}/ParametricMap.cpp"
    )
    target_include_directories(
        testAdvection
        PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}"
    )
    target_link_libraries(testAdvection LINK_PUBLIC doctest::doctest Eigen3::Eigen)

    add_executable(
        testAdvectionPeriodicBC
        AdvectionPeriodicBC_test.cpp
        "${SRC_DIR}/DGTransport.cpp"
        "${SRC_DIR}/ParametricMap.cpp"
        "${SRC_DIR}/Interpolations.cpp"
    )
    target_include_directories(
        testAdvectionPeriodicBC
        PRIVATE "${CoreDir}" "${SRC_DIR}" "${CoreDir}/${ModelArrayStructure}"
    )
    target_link_libraries(testAdvectionPeriodicBC LINK_PUBLIC doctest::doctest Eigen3::Eigen)

    set(TEST_EXECUTABLES
        testDGModelArray
        testDGVectorHolder
        testCGModelArray
        testParametricMesh
        testParametricMeshArea
        testAdvection
        testAdvectionPeriodicBC
    )

    foreach(test_exec IN LISTS TEST_EXECUTABLES)
        add_test(NAME ${test_exec} COMMAND ${test_exec})
    endforeach()
endif()
