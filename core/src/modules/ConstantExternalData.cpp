/*!
 * @file ConstantExternalData.cpp
 *
 * @date Mar 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ConstantExternalData.hpp"

namespace Nextsim {

const double ConstantExternalData::c_t_air = -1;
const double ConstantExternalData::c_dew_2m = -4;
const double ConstantExternalData::c_p_air = 1e5;
const double ConstantExternalData::c_mixrat = -1; // invalid: use dew point instead
const double ConstantExternalData::c_swin = 0; // night
const double ConstantExternalData::c_lwin = 311;
const double ConstantExternalData::c_mld = 10;
const double ConstantExternalData::c_snowfall = 0;

void ConstantExternalData::getData(const TimePoint& tp)
{
    // Make sure the arrays have the correct size
    t_air.resize();
    dew_2m.resize();
    p_air.resize();
    mixrat.resize();
    swin.resize();
    lwin.resize();
    mld.resize();
    snowfall.resize();

    // Fill the arrays with data. Yes, they are probably all already full of
    // this data, but if you are using this module, is performance really that
    // important to you?
    // H arrays
    for (size_t i = 0; i < t_air.size(); ++i) {
        t_air[i] = c_t_air;
        dew_2m[i] = c_dew_2m;
        c_p_air[i] = c_p_air;
        mixrat[i] = c_mixrat;
        swin[i] = c_swin;
        lwin[i] = c_lwin;
        mld[i] = c_mld;
        snowfall[i] = c_snowfall;
    }
}

} /* namespace Nextsim */
