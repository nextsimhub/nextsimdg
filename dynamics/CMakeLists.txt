add_subdirectory(src)

if (BUILD_TESTS)
    add_subdirectory(test)
endif()

set(NextsimSources
    "${NextsimSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    PARENT_SCOPE)
