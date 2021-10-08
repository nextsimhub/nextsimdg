# Build the unit, integration and model tests for neXtSIM

set(SRC_DIR "../src")
set(INCLUDE_DIR "${SRC_DIR}/include")

# add_executable(testexe
#   test/TestSrc.cpp
#   otherSource.cpp)
#target_link_libraries(testexe PRIVATE Catch2::Catch2)

# include(CTest)
# include(Catch)
# catch_discover_tests(testexe)

include_directories(${INCLUDE_DIR})

add_executable(testTimed
    "Timed_test.cpp"
    "${SRC_DIR}/CountTimer.cpp"
    "${SRC_DIR}/Timed.cpp"
    "${INCLUDE_DIR}/Timed.hpp"
    "${INCLUDE_DIR}/Timer.hpp"
    )
target_link_libraries(testTimed PRIVATE Catch2::Catch2)

add_executable(testIterator
    "Iterator_test.cpp"
    "${SRC_DIR}/Iterator.cpp"
    "${SRC_DIR}/Timed.cpp"
    "${SRC_DIR}/CountTimer.cpp"
    "${SRC_DIR}/Logged.cpp"
    "${INCLUDE_DIR}/Iterator.hpp"
    "${INCLUDE_DIR}/Timed.hpp"
    "${INCLUDE_DIR}/Logged.hpp"
    )
target_link_libraries(testIterator PRIVATE Catch2::Catch2)

add_executable(testSimpleIter
    "SimpleIterant_test.cpp"
    "${SRC_DIR}/SimpleIterant.cpp"
    "${SRC_DIR}/Iterator.cpp"
    "${SRC_DIR}/Timed.cpp"
    "${SRC_DIR}/CountTimer.cpp"
    "${SRC_DIR}/Logged.cpp"
    "${INCLUDE_DIR}/SimpleIterant.hpp"
    "${INCLUDE_DIR}/Iterator.hpp"
    "${INCLUDE_DIR}/Timed.hpp"
    "${INCLUDE_DIR}/Logged.hpp"
    )
target_link_libraries(testSimpleIter PRIVATE Catch2::Catch2)

add_executable(testCommandLineParser
    "CommandLineParser_test.cpp"
    "${SRC_DIR}/CommandLineParser.cpp"
    )
target_include_directories(testCommandLineParser PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(testCommandLineParser LINK_PUBLIC ${Boost_LIBRARIES} Catch2::Catch2)
