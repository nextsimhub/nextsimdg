/*
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Örn Ólason <einar.olason@nersc.no>
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
    IIceAlbedo()
        : iceAlbedoAccessor(getStore(), RO, ModelArray::Type::H)
        , icePenSWAccessor(getStore(), RO, ModelArray::Type::H)
        , hsnowAccessor(getStore())
        , ciceAccessor(getStore())
    {
    }
    virtual ~IIceAlbedo() = default;

    std::string getName () const { return "IIceAlbedo"; }
    
    /*!
     * @brief A virtual function that calculates the ice surface short-wave
     * albedo and fraction of penetrating short-wave radiation.
     *
     * @param temperature The temperature of the ice surface.
     * @param snowThickness The true snow thickness on top of the ice.
     * @param i0 The fraction of short-wave radiation that can penetrate bare ice (not taking snow
     * cover into account).
     */
    virtual std::tuple<double, double> surfaceShortWaveBalance(
        double temperature, double snowThickness, double i0)
        = 0;

    virtual void setData(const ModelState::DataMap&)
    {
        iceAlbedoAccessor.getHostRW().resize();
        icePenSWAccessor.getHostRW().resize();
    }

    virtual void update(const TimestepTime& ts) { }

    /*!
     * Sets the time parameter for the implementation, if it is time dependent.
     * @param time The desired TimePoint.
     */
    virtual void setTime(const TimePoint& tp) { }

protected:
    // owned fields
    ModelArrayAccessor<Protected::ICE_ALBEDO, RW> iceAlbedoAccessor;
    ModelArrayAccessor<Protected::ICE_PEN_SW, RW> icePenSWAccessor;

    // needed for snow thickness
    ModelArrayAccessor<Shared::H_SNOW_DG> hsnowAccessor;
    ModelArrayAccessor<Shared::C_ICE_DG> ciceAccessor;
};
}
#endif /* IICEALBEDO_HPP */
