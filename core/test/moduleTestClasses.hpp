/*
 * @file moduleTestClasses.hpp
 *
 * @date Oct 21, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef TEST_TESTCLASSES_HPP
#define TEST_TESTCLASSES_HPP

// Test classes
class ITest {
public:
    virtual ~ITest() = default;
    virtual int operator()() { return 0; }
};
class Impl1 : public ITest {
public:
    ~Impl1() = default;
    int operator()() override { return 1; }
};
class Impl2 : public ITest {
public:
    ~Impl2() = default;
    int operator()() override { return 2; }
};

#endif /* TEST_TESTCLASSES_HPP */
