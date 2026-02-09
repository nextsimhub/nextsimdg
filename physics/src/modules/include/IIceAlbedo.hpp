/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Örn Ólason <einar.olason@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef IICEALBEDO_HPP
#define IICEALBEDO_HPP

#include "include/ModelArrayAccessor.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"
#include <tuple>

namespace Nextsim {
//! The interface class for ice albedo calculation.
class IIceAlbedo : public ModelComponent {
public:
    IIceAlbedo();

    virtual ~IIceAlbedo() = default;

    virtual std::string getName() const = 0;

    virtual void setData(const ModelState::DataMap&);

    /*!
     * @brief Calculates the ice surface short-wave albedo and fraction of penetrating short-wave
     * radiation and writes the results to the ICE_ALBEDO and ICE_PEN_SW fields.
     */
    virtual void update(const TimestepTime& ts) = 0;

protected:
    // owned fields
    ModelArrayAccessor<Protected::ICE_ALBEDO, RW>
        iceAlbedoAccessor; // ice surface short-wave albedo
    ModelArrayAccessor<Protected::ICE_PEN_SW, RW>
        icePenSWAccessor; // fraction of penetrating short-wave radiation

    // needed for snow thickness
    ModelArrayAccessor<Shared::H_SNOW_DG> hsnowAccessor;
    ModelArrayAccessor<Shared::C_ICE_DG> ciceAccessor;
    ModelArrayAccessor<Protected::T_SURF> tsurfAccessor; // temperature of the ice surface
};
}
#endif /* IICEALBEDO_HPP */
