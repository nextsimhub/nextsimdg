set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

add_subdirectory("${ModuleDir}")

set(Sources
    "BenchmarkCoordinates.cpp"
    "IceGrowth.cpp"
    "IceMinima.cpp"
    "NetCDFForcings.cpp"
    "SlabOcean.cpp"
    )

list(TRANSFORM Sources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(NextsimSources "${NextsimSources}" "${Sources}" PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE
)
