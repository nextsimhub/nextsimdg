add_subdirectory(src)
#if (NOT(CMAKE_BUILD_TYPE STREQUAL "Release"))
#    add_subdirectory(test)
#endif()

add_subdirectory(applications)
add_subdirectory(tests)


set(NextsimSources
    "${NextsimSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    PARENT_SCOPE)