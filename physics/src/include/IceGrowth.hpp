/*!
 * @file IceGrowth.hpp
 *
 * @date 04 Jun 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef ICEGROWTH_HPP
#define ICEGROWTH_HPP

#include "include/Configured.hpp"
#include "include/IDamageHealing.hpp"
#include "include/IIceThermodynamics.hpp"
#include "include/ILateralIceSpread.hpp"
#include "include/IceMinima.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class IceGrowth : public ModelComponent, public Configured<IceGrowth> {
public:
    IceGrowth();
    virtual ~IceGrowth() = default;

    enum {
        ICE_THERMODYNAMICS_KEY,
        LATERAL_GROWTH_KEY,
        MINC_KEY,
        MINH_KEY,
        USE_THERMO_KEY,
    };

    void configure() override;
    ConfigMap getConfiguration() const override;

    std::string getName() const override { return "IceGrowth"; }

    void setData(const ModelState::DataMap&) override;
    ModelState getStateDiagnostic() const override;
    ModelState getStatePrognostic() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void update(const TimestepTime&);

    /*!
     * Updates the true ice and snow thickness arrays from the cell averages.
     */
    void initializeThicknesses();

private:
    // Vertical Growth ModelComponent & Module
    std::unique_ptr<IIceThermodynamics> iVertical;
    // Lateral Growth ModuleComponent & Module
    std::unique_ptr<ILateralIceSpread> iLateral;
    // Damage Healing ModuleComponent & Module
    std::unique_ptr<IDamageHealing> iHealing;

    // Data fields
    // Owned, shared data fields
    HField hice; // Updated true ice thickness, m
    HField cice; // Updated ice concentration
    HField hsnow; // Updated true snow thickness, m
    HField hice0; // Timestep initial true ice thickness, m
    HField hsnow0; // Timestep initial true snow thickness, m

    ModelArrayRef<Protected::H_ICE> hIceCell; // Timestep initial cell averaged ice thickness, m
    ModelArrayRef<Protected::H_SNOW> hSnowCell; // Timestep initial cell averaged snow thickness, m
    ModelArrayRef<Protected::C_ICE> cice0; // Timestep initial ice concentration
    ModelArrayRef<Shared::Q_OW, RW> qow; // open water heat flux, from FluxCalculation
    ModelArrayRef<Shared::DELTA_HICE> deltaHi; // New ice thickness this timestep, m

    bool doThermo = true; // Perform any thermodynamics calculations at all

    void initializeThicknessesElement(size_t i, const TimestepTime&);
};

} /* namespace Nextsim */

#endif /* ICEGROWTH_HPP */
