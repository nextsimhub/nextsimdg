# Include dynamics headers
set(DynamicsSrc "../src")
set(SRC_DIR "${DynamicsSrc}")
set(CoreSrc "../../core/src")
set(Core_DIR "${CoreSrc}")
set(INCLUDE_DIR "${DynamicsSrc}/include" "${CoreSrc}/include")

set(NEXTSIMLIB nextsimlib)

#include_directories(${INCLUDE_DIR})

add_executable(benchmark_mehlmann_sasip_mevp benchmark_mehlmann_sasip_mevp.cpp
          "${SRC_DIR}/stopwatch.cpp"
          "${SRC_DIR}/DGTransport.cpp"
          "${SRC_DIR}/Interpolations.cpp"
          "${SRC_DIR}/ParametricMesh.cpp"
          "${SRC_DIR}/ParametricMap.cpp"
          "${SRC_DIR}/cgParametricMomentum.cpp"
          "${SRC_DIR}/VectorManipulations.cpp"
          )
target_link_libraries(benchmark_mehlmann_sasip_mevp LINK_PUBLIC
          ${NEXTSIMLIB}
	  Eigen3::Eigen)

target_include_directories(benchmark_mehlmann_sasip_mevp
    PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${NextsimIncludeDirs}"
)
