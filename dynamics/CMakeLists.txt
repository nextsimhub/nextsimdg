add_subdirectory(src)
#if (NOT(CMAKE_BUILD_TYPE STREQUAL "Release"))
#    add_subdirectory(test)
#endif()

add_subdirectory(applications)
add_subdirectory(tests)

