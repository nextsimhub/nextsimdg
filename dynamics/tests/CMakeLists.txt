# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")
set(CoreDir "../../core/src/")

include_directories(${INCLUDE_DIR})

add_executable(example1-parametricmesh example1-parametricmesh.cpp
"${SRC_DIR}/stopwatch.cpp"
"${SRC_DIR}/parametricTransport.cpp"
"${SRC_DIR}/Interpolations.cpp"
"${SRC_DIR}/ParametricMesh.cpp"
"${SRC_DIR}/parametricTools.cpp" 
"${SRC_DIR}/cgParametricMomentum.cpp" 
)
target_include_directories(example1-parametricmesh PRIVATE "${CoreDir}")
target_link_libraries( example1-parametricmesh LINK_PUBLIC Eigen3::Eigen)

add_executable(dgma_test
    "DGModelArray_test.cpp"
    "${CoreDir}/ModelArray.cpp"
    )
target_include_directories(dgma_test PRIVATE "${CoreDir}" "${SRC_DIR}")
target_link_libraries(dgma_test LINK_PUBLIC Catch2::Catch2 Eigen3::Eigen)

add_executable(cgma_test
    "CGModelArray_test.cpp"
    "${CoreDir}/ModelArray.cpp"
    )
target_include_directories(cgma_test PRIVATE "${CoreDir}" "${SRC_DIR}")
target_link_libraries(cgma_test LINK_PUBLIC Catch2::Catch2 Eigen3::Eigen)