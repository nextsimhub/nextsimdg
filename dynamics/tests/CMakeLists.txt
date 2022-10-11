# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")

include_directories(${INCLUDE_DIR})

add_executable(example1-parametricmesh example1-parametricmesh.cpp
"${SRC_DIR}/stopwatch.cpp"
"${SRC_DIR}/parametricTransport.cpp"
"${SRC_DIR}/Interpolations.cpp"
"${SRC_DIR}/ParametricMesh.cpp"
"${SRC_DIR}/parametricTools.cpp" 
"${SRC_DIR}/cgParametricMomentum.cpp" 
)

target_link_libraries( example1-parametricmesh LINK_PUBLIC Eigen3::Eigen)
