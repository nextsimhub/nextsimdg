/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MISSINGDATA_HPP
#define MISSINGDATA_HPP

#include "include/Configured.hpp"
#include "include/FloatType.hpp"

namespace Nextsim {

class MissingData : public Configured<MissingData> {
public:
    static constexpr FloatType defaultValue = 1.7e38;
    inline static FloatType value() { return getValue(); }
    inline static void setValue(FloatType mdi) { getValue() = mdi; }

private:
    inline static FloatType& getValue()
    {
        static FloatType value = defaultValue;
        return value;
    }
};

} /* namespace Nextsim */

#endif /* MISSINGDATA_HPP */
