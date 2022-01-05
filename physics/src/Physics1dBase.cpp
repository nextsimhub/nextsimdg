/*!
 * @file Physics1dBase.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Physics1dBase.hpp"

namespace Nextsim {

void Physics1dBase::physics1d(ElementData& data)
{
    data.updateDerivedData(data, data, data);
    data.calculate(data, data, data);
}
} /* namespace Nextsim */
