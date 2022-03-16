/*!
 * @file IExternalData.hpp
 *
 * @date Mar 10, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/Time.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelModule.hpp"

#ifndef IEXTERNALDATA_HPP
#define IEXTERNALDATA_HPP

namespace Nextsim {
class IExternalData : public ModelModule {
public:
    IExternalData()
    {
        registerProtectedArray(ProtectedArray::T_AIR, &t_air);
        registerProtectedArray(ProtectedArray::DEW_2M, &dew_2m);
        registerProtectedArray(ProtectedArray::P_AIR, &p_air);
        registerProtectedArray(ProtectedArray::MIXRAT, &mixrat);
        registerProtectedArray(ProtectedArray::SW_IN, &swin);
        registerProtectedArray(ProtectedArray::LW_IN, &lwin);
        registerProtectedArray(ProtectedArray::MLD, &mld);
        registerProtectedArray(ProtectedArray::SNOW, &snowfall);
    }
    virtual ~IExternalData() = default;

    virtual void getData(const TimePoint&);
private:
    HField t_air;
    HField dew_2m;
    HField p_air;
    HField mixrat;
    HField swin;
    HField lwin;
    HField mld;
    HField snowfall;
};

}

#endif /* IEXTERNALDATA_HPP */
