/*!
 * @file DummyExternalData.hpp
 *
 * @date Jan 20, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_DUMMYEXTERNALDATA_HPP
#define CORE_SRC_INCLUDE_DUMMYEXTERNALDATA_HPP

#include "include/ExternalData.hpp"
#include "include/IStructure.hpp"

namespace Nextsim {

//! A class to assign fixed, constant values to the ExternalData members of an
//! IStructure.
class DummyExternalData {
public:
    ~DummyExternalData() = default;

    static void setAll(IStructure& is)
    {
        for (is.cursor = 0; is.cursor; ++is.cursor) {
            is.cursor->airTemperature() = -1;
            is.cursor->dewPoint2m() = -4;
            is.cursor->airPressure() = 1e5;
            is.cursor->ExternalData::mixingRatio() = -1.;
            is.cursor->incomingShortwave() = 0; // night
            is.cursor->incomingLongwave() = 311;
            is.cursor->mixedLayerDepth() = 10;
            is.cursor->snowfall() = 0;
        }
    }

private:
    DummyExternalData() = default;
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DUMMYEXTERNALDATA_HPP */
