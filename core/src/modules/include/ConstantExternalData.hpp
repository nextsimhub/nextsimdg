/*!
 * @file ConstantExternalData.hpp
 *
 * @date Mar 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_MODULES_INCLUDE_CONSTANTEXTERNALDATA_HPP
#define CORE_SRC_MODULES_INCLUDE_CONSTANTEXTERNALDATA_HPP

#include "include/IExternalData.hpp"

namespace Nextsim {

class ConstantExternalData : public IExternalData {
public:
    ConstantExternalData() = default;
    virtual ~ConstantExternalData() = default;

    void getData(const TimePoint&) override;
private:
    // Hardcoded values for the external data
    static const double c_t_air;
    static const double c_dew_2m;
    static const double c_p_air;
    static const double c_mixrat;
    static const double c_swin;
    static const double c_lwin;
    static const double c_mld;
    static const double c_snowfall;

};

} /* namespace Nextsim */

#endif /* CORE_SRC_MODULES_INCLUDE_CONSTANTEXTERNALDATA_HPP */
