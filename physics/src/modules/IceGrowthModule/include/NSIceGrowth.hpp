/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#ifndef NSICEGROWTH_HPP
#define NSICEGROWTH_HPP

#include "include/IIceGrowth.hpp"

namespace Nextsim {

class NSIceGrowth : public IIceGrowth, public Configured<NSIceGrowth> {
public:
    NSIceGrowth();
    virtual ~NSIceGrowth() = default;

    enum {
        ICE_THERMODYNAMICS_KEY,
        LATERAL_GROWTH_KEY,
        MINC_KEY,
        MINH_KEY,
        USE_THERMO_KEY,
    };

    void configure() override;
    ConfigMap getConfiguration() const override;

    std::string getName() const override { return "NSIceGrowth"; }

    void setData(const ModelState::DataMap&) override;
    ModelState getStateDiagnostic() const override;
    ModelState getStatePrognostic() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void update(const TimestepTime&) override;

private:
    // Vertical Growth ModelComponent & Module
    std::unique_ptr<IIceThermodynamics> iVertical;
    // Lateral Growth ModuleComponent & Module
    std::unique_ptr<ILateralIceSpread> iLateral;
    // Damage Healing ModuleComponent & Module
    std::unique_ptr<IDamageHealing> iHealing;

    // Data fields
    ModelArrayRef<Shared::H_ICE_DG> hice; // Timestep initial cell averaged ice thickness, m
    ModelArrayRef<Shared::H_SNOW_DG> hsnow; // Timestep initial cell averaged snow thickness, m
    ModelArrayRef<Shared::C_ICE_DG> cice; // Timestep initial ice concentration
    ModelArrayRef<Shared::Q_OW, RW> qow; // open water heat flux, from FluxCalculation
    ModelArrayRef<Shared::DELTA_HICE> deltaHi; // New ice thickness this timestep, m

    bool doThermo = true; // Perform any thermodynamics calculations at all
};

} /* namespace Nextsim */

#endif /* NSICEGROWTH_HPP */
