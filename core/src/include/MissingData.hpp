/*!
 * @file MissingData.hpp
 *
 * @date Jun 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MISSINGDATA_HPP
#define MISSINGDATA_HPP

#include "include/Configured.hpp"

namespace Nextsim {

class MissingData : public Configured<MissingData> {
public:
    static constexpr double defaultValue = 1.7e38;
    inline static double value() { return getValue(); }
    inline static void setValue(double mdi) { getValue() = mdi; }

private:
    inline static double& getValue()
    {
        static double value = defaultValue;
        return value;
    }
};

} /* namespace Nextsim */

#endif /* MISSINGDATA_HPP */
