# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include")

include_directories(${INCLUDE_DIR})

#   add_executable(example1-sasipmesh example1-sasipmesh.cpp
#   "${SRC_DIR}/stopwatch.cpp"
#   "${SRC_DIR}/parametricTransport.cpp"
#   "${SRC_DIR}/Interpolations.cpp"
#   "${SRC_DIR}/ParametricMesh.cpp"
#   "${SRC_DIR}/parametricTools.cpp" 
#   "${SRC_DIR}/cgParametricMomentum.cpp" 
#   )
#   
#   
#   add_executable(example2-sasipmesh example2-sasipmesh.cpp
#   "${SRC_DIR}/stopwatch.cpp"
#   "${SRC_DIR}/parametricTransport.cpp"
#   "${SRC_DIR}/Interpolations.cpp"
#   "${SRC_DIR}/ParametricMesh.cpp"
#   "${SRC_DIR}/parametricTools.cpp" 
#   "${SRC_DIR}/cgParametricMomentum.cpp" 
#   )
#   
#   add_executable(example3-orcamesh example3-orcamesh.cpp
#   "${SRC_DIR}/stopwatch.cpp"
#   "${SRC_DIR}/parametricTransport.cpp"
#   "${SRC_DIR}/Interpolations.cpp"
#   "${SRC_DIR}/ParametricMesh.cpp"
#   "${SRC_DIR}/parametricTools.cpp" 
#   "${SRC_DIR}/cgParametricMomentum.cpp" 
#   )
#   

add_executable(example4-orcamesh-spherical example4-orcamesh-spherical.cpp
"${SRC_DIR}/stopwatch.cpp"
"${SRC_DIR}/SphericalTransport.cpp"
"${SRC_DIR}/Interpolations.cpp"
"${SRC_DIR}/ParametricMesh.cpp"
"${SRC_DIR}/parametricTools.cpp"
"${SRC_DIR}/SphericalTools.cpp" 
"${SRC_DIR}/cgParametricMomentum.cpp" 
)

add_executable(example4-orcamesh-spherical-greenland example4-orcamesh-spherical-greenland.cpp
"${SRC_DIR}/stopwatch.cpp"
"${SRC_DIR}/SphericalTransport.cpp"
"${SRC_DIR}/Interpolations.cpp"
"${SRC_DIR}/ParametricMesh.cpp"
"${SRC_DIR}/parametricTools.cpp"
"${SRC_DIR}/SphericalTools.cpp" 
"${SRC_DIR}/cgParametricMomentum.cpp" 
)

#  target_link_libraries(example1-sasipmesh LINK_PUBLIC Eigen3::Eigen)
#  target_link_libraries(example2-sasipmesh LINK_PUBLIC Eigen3::Eigen)
#  target_link_libraries(example3-orcamesh  LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(example4-orcamesh-spherical  LINK_PUBLIC Eigen3::Eigen)
target_link_libraries(example4-orcamesh-spherical-greenland  LINK_PUBLIC Eigen3::Eigen)
