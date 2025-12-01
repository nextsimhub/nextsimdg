/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#ifndef IICEGROWTH_HPP
#define IICEGROWTH_HPP

#include "include/Configured.hpp"
#include "include/IDamageHealing.hpp"
#include "include/IIceThermodynamics.hpp"
#include "include/ILateralIceSpread.hpp"
#include "include/IceMinima.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class IIceGrowth : public ModelComponent {
public:
    virtual ~IIceGrowth() = default;

    virtual void update(const TimestepTime&) = 0;
protected:
    IIceGrowth() = default;

    static void setCMin(double cMin) { IceMinima::cMin = cMin; };
    static void setHMin(double hMin) { IceMinima::hMin = hMin; };
};

} /* namespace Nextsim */

#endif /* IICEGROWTH_HPP */
