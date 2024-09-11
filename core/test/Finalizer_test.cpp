/*!
 * @file ModelFinalizer_test.cpp
 *
 * @date Sep 3, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/Finalizer.hpp"

namespace Nextsim {

// A mock Model class. Does nothing except call copy the behaviour of the real
// Model class in its destructor
class Model
{
public:
    Model() = default;
    ~Model() { Finalizer::finalize(); }
};

class TheCount
{
public:
    inline static size_t increment() { return get()++; }
    inline static size_t preincrement() { return ++get(); }
    inline static void incrementAndPrint() { std::cout << "TheCount=" << preincrement() << std::endl; }
    inline static size_t value() { return get(); }
    inline static size_t count() { return get(); }
    inline static void reset() { get() = 0; }
protected:
    inline static size_t& get()
    {
        static size_t count;
        return count;
    }
};

// In the British aristocracy an Earl is equivalent to a count in the rest of Europe!
class Earl
{
public:
    inline static void incr() { TheCount::increment(); }
};

class HappyTime : public std::exception
{
public:
    const char* what() const noexcept { return "Nothing is wrong :D"; }
};

void throwHappy() { throw HappyTime(); }

TEST_SUITE_BEGIN("ModelFinalizer");
TEST_CASE("Duplicate functions")
{
    TheCount::reset();
    REQUIRE(TheCount::count() == 0);
//    Finalizer::clear();
    REQUIRE(Finalizer::count() == 0);
    Finalizer::registerFunc(TheCount::increment);
    Finalizer::registerFunc(TheCount::increment);
    REQUIRE(Finalizer::count() == 2);
    Finalizer::finalize();
    REQUIRE(TheCount::count() == 2);
}

TEST_CASE("Unique functions")
{
    TheCount::reset();
    REQUIRE(TheCount::count() == 0);
//    Finalizer::clear();
    REQUIRE(Finalizer::count() == 0);
    Finalizer::registerUnique(TheCount::increment);
    REQUIRE(Finalizer::contains(TheCount::increment));
    Finalizer::registerUnique(TheCount::increment);
    REQUIRE(Finalizer::count() == 1);
    Finalizer::finalize();
    REQUIRE(TheCount::count() == 1);
}

TEST_CASE("Finalize")
{
    TheCount::reset();
    REQUIRE(TheCount::count() == 0);
//    Finalizer::clear();
    REQUIRE(Finalizer::count() == 0);
    Finalizer::registerFunc(TheCount::increment);
    Finalizer::registerFunc(TheCount::increment);
    REQUIRE(Finalizer::count() == 2);
    Finalizer::finalize();
    REQUIRE(Finalizer::count() == 0);
    REQUIRE(TheCount::count() == 2);
    TheCount::reset();
    Finalizer::finalize();
    REQUIRE(Finalizer::count() == 0);
    REQUIRE(TheCount::count() == 0);
}

TEST_CASE("Exception")
{
    TheCount::reset();
    REQUIRE(TheCount::count() == 0);
//    Finalizer::clear();
    REQUIRE(Finalizer::count() == 0);
    Finalizer::registerFunc(TheCount::increment);
    Finalizer::registerFunc(throwHappy);
    Finalizer::registerFunc(TheCount::increment);
    REQUIRE(Finalizer::count() == 3);
    REQUIRE_THROWS_AS(Finalizer::finalize(), HappyTime);
    // Only one increment occurred because of the exception
    REQUIRE(TheCount::count() == 1);
    // The executed increment and throwing function have been removed. Only one function (the
    // second increment) should remain.
    REQUIRE(Finalizer::count() == 1);
}
} /* namespace Nextsim */
