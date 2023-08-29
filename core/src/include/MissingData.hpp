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
    static const double defaultValue;
    static double value;
};

} /* namespace Nextsim */

#endif /* MISSINGDATA_HPP */
