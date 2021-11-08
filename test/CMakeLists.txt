# Build the unit, integration and model tests for neXtSIM

set(SRC_DIR "../src")
set(INCLUDE_DIR "${SRC_DIR}/include")

# add_executable(testexe
#	test/TestSrc.cpp
#	otherSource.cpp)
#target_link_libraries(testexe PRIVATE Catch2::Catch2)

# include(CTest)
# include(Catch)
# catch_discover_tests(testexe)

include_directories(${INCLUDE_DIR})

add_executable(testIterator
    "Iterator_test.cpp"
    "${SRC_DIR}/Iterator.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/Logged.cpp"
    )
target_link_libraries(testIterator PRIVATE Catch2::Catch2)

add_executable(testSimpleIter
    "SimpleIterant_test.cpp"
    "${SRC_DIR}/SimpleIterant.cpp"
    "${SRC_DIR}/Iterator.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/Logged.cpp"
    )
target_link_libraries(testSimpleIter PRIVATE Catch2::Catch2)

add_executable(testTimer
    "Timer_test.cpp"
    "${SRC_DIR}/Timer.cpp"
    )
target_link_libraries(testTimer PRIVATE Catch2::Catch2)
target_include_directories(testTimer PRIVATE "${SRC_DIR}")

add_executable(testScopedTimer
    "ScopedTimer_test.cpp"
    "${SRC_DIR}/Timer.cpp"
    "${SRC_DIR}/ScopedTimer.cpp"
    )
target_link_libraries(testScopedTimer PRIVATE Catch2::Catch2)
target_include_directories(testScopedTimer PRIVATE "${SRC_DIR}")


