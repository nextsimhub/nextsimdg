/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#ifndef NSCOLUMNPHYSICS_HPP
#define NSCOLUMNPHYSICS_HPP

#include "include/IColumnPhysics.hpp"

namespace Nextsim {

/*!
 * The column physics package derived from neXtSIM-Lagrangian
 */
class NSColumnPhysics : public IColumnPhysics, public Configured<NSColumnPhysics> {
public:
    NSColumnPhysics();
    virtual ~NSColumnPhysics() = default;

    enum {
        ICE_THERMODYNAMICS_KEY,
        LATERAL_GROWTH_KEY,
        MINC_KEY,
        MINH_KEY,
        USE_THERMO_KEY,
    };

    void configure() override;
    ConfigMap getConfiguration() const override;

    std::string getName() const override { return "NSColumnPhysics"; }

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
    ModelArrayAccessor<Shared::H_ICE_DG>
        hiceAccessor; // Timestep initial cell averaged ice thickness, m
    ModelArrayAccessor<Shared::H_SNOW_DG>
        hsnowAccessor; // Timestep initial cell averaged snow thickness, m
    ModelArrayAccessor<Shared::C_ICE_DG> ciceAccessor; // Timestep initial ice concentration
    ModelArrayAccessor<Shared::Q_OW, RW> qowAccessor; // open water heat flux, from FluxCalculation
    ModelArrayAccessor<Shared::DELTA_HICE> deltaHiAccessor; // New ice thickness this timestep, m

    bool doThermo = true; // Perform any thermodynamics calculations at all
};

} /* namespace Nextsim */

#endif /* NSCOLUMNPHYSICS_HPP */
