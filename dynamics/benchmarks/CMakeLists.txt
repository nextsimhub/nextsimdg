# Standalone executables of the Mehlmann et al. (2021) benchmark, running the dynamics kernels
# without the rest of the model. Built with -DBUILD_BENCHMARKS=ON.

set(BenchmarkNames benchmark_mehlmann_mevp benchmark_mehlmann_bbm benchmark_mehlmann_ldg_mevp)

foreach(name IN LISTS BenchmarkNames)
    add_executable(${name} "${name}.cpp")
    target_link_libraries(${name} PRIVATE nextsimlib)
    target_include_directories(
        ${name} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}" ${NextsimIncludeDirs})
    target_compile_definitions(${name} PRIVATE DGCOMP=${DGComp} DGSTRESSCOMP=${DGStressComp}
                                               CGDEGREE=${CGDegree} LDGCOMP=${LDGComp})
endforeach()
