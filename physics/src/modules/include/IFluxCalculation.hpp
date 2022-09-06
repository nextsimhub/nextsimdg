/*!
 * @file IFluxCalculation.hpp
 *
 * @date Apr 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IFLUXCALCULATION_HPP
#define IFLUXCALCULATION_HPP

#include "include/AtmosphereOceanState.hpp"
#include "include/Configured.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"

namespace Nextsim {
//! An interface class for calculating ocean-ice-atmosphere fluxes
class IFluxCalculation : public ModelComponent, public Configured<IFluxCalculation> {
public:
    IFluxCalculation()
        : qow(ModelArray::Type::H)
        , subl(ModelArray::Type::H)
        , qia(ModelArray::Type::H)
        , dqia_dt(ModelArray::Type::H)
        , qio(ModelArray::Type::H)
    {
        // register shared arrays
        registerSharedArray(SharedArray::Q_OW, &qow);
        registerSharedArray(SharedArray::SUBLIM, &subl);
        registerSharedArray(SharedArray::Q_IA, &qia);
        registerSharedArray(SharedArray::DQIA_DT, &dqia_dt);
        registerSharedArray(SharedArray::Q_IO, &qio);
    }
    virtual ~IFluxCalculation() = default;

    void setData(const ModelState::DataMap& ms) override
    {
        aoState.setData(ms);

        qow.resize();
        subl.resize();
        qia.resize();
        dqia_dt.resize();
        qio.resize();
    }

    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        std::cout << "IFluxCalculation::getStateRecursive: data.size=" << 0 << ", config.size=" << 0
                  << std::endl;
        return os ? getState() : ModelState();
    }

    std::string getName() const override { return "IFluxCalculation"; }

    std::unordered_set<std::string> hFields() const override
    {
        return { "qow", "subl", "qia", "dqia_dt", "qio" };
    }

    void configure() override { aoState.configure(); }
    /*!
     * Updates the flux calculation for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    virtual void update(const TimestepTime&) = 0;

protected:
    // All fluxes are positive upwards, including incident radiation fluxes
    // Owned, shared fields
    HField qow; // Open water heat flux [W m⁻²]
    HField subl; // Ice sublimative mass flux [kg m⁻²]
    HField qia; // Ice-atmosphere heat flux [W m⁻²]
    HField dqia_dt; // Derivative of qia w.r.t. ice surface temperature
    HField qio; // Ice-ocean heat flux [W m⁻²]

    // The interface to the rest of the world
    AtmosphereOceanState aoState;
};
}
#endif /* IFLUXCALCULATION_HPP */
