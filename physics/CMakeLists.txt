add_subdirectory(src)
if (NOT(CMAKE_BUILD_TYPE STREQUAL "Release"))
    add_subdirectory(test)
endif()

set(NextsimSources
    "${NextsimSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    PARENT_SCOPE)
