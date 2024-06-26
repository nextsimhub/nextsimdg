/*!
 * @file OutputSpec.hpp
 *
 * @date Aug 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef OUTPUTSPEC_HPP
#define OUTPUTSPEC_HPP

#include <bitset>

namespace Nextsim {
// TODO: Replace this with a more fine grained output specification class
class OutputSpec {
public:
    OutputSpec()
    {
        bitVector.reset();
    }
    OutputSpec(bool activate)
        : OutputSpec()
    {
        if (!activate) bitVector.set(SUPPRESS_OUTPUT);
    }

    void setSuppressOutput() { bitVector.set(SUPPRESS_OUTPUT); }
    bool suppressOutput() const { return bitVector.test(SUPPRESS_OUTPUT); }
    operator bool() const { return !suppressOutput(); }
    void setAllComponents() { bitVector.set(ALL_COMPONENTS); }
    bool allComponents() const { return bitVector.test(ALL_COMPONENTS); }
private:
    enum {
        SUPPRESS_OUTPUT,
        ALL_COMPONENTS,
        OUTPUTSPEC_COUNT
    };

    std::bitset<OUTPUTSPEC_COUNT> bitVector;
};
}

#endif /* OUTPUTSPEC_HPP */
