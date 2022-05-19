/*!
 * @file IceGrowth.hpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PHYSICS_SRC_INCLUDE_ICEGROWTH_HPP
#define PHYSICS_SRC_INCLUDE_ICEGROWTH_HPP

#include "include/Configured.hpp"
#include "include/IFluxCalculation.hpp"
#include "include/IIceThermodynamics.hpp"
#include "include/ILateralIceSpread.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class IceGrowth : public ModelComponent, public Configured<IceGrowth> {
public:
    IceGrowth();
    virtual ~IceGrowth() = default;

    void configure() override;
    enum {
        ICE_THERMODYNAMICS_KEY,
        LATERAL_GROWTH_KEY,
        FLUX_CALCULATOR_KEY,
        MINC_KEY,
        MINH_KEY,
    };

    std::string getName() const override { return "IceGrowth"; }

    void setData(const ModelState&) override;
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::set<std::string> hFields() const override
    {
        return { "updated_hice", "updated_cice", "updated_hsnow" };
    }
    std::set<std::string> uFields() const override { return {}; }
    std::set<std::string> vFields() const override { return {}; }
    std::set<std::string> zFields() const override { return {}; }

    void update(const TimestepTime&);

    static double minimumIceThickness() { return minh; }
    static double minimumIceConcentration() { return minc; }

private:
    // Vertical Growth ModelComponent & Module
    std::unique_ptr<IIceThermodynamics> iVertical;
    // Lateral Growth ModuleComponent & Module
    std::unique_ptr<ILateralIceSpread> iLateral;
    // Flux calculation Module Component
    std::unique_ptr<IFluxCalculation> iFluxes;

    // Data fields
    // Owned, shared data fields
    HField hice; // Updated true ice thickness, m
    HField cice; // Updated ice concentration
    HField hsnow; // Updated true snow thickness, m
    HField newice; // New ice over open water this timestep, m

    // Owned data fields, not shared
    HField deltaCFreeze; // New ice concentration due to freezing (+ve)
    HField deltaCMelt; // Ice concentration loss due to melting (-ve)

    ModelArrayRef<ProtectedArray::HTRUE_ICE> hice0; // initial ice thickness (ice averaged)
    ModelArrayRef<ProtectedArray::C_ICE> cice0; // initial ice concentration
    ModelArrayRef<ProtectedArray::HTRUE_SNOW> hsnow0; // initial snow thickness (ice averaged)
    ModelArrayRef<SharedArray::Q_OW, RW> qow; // open water heat flux, from FluxCalculation
    ModelArrayRef<ProtectedArray::ML_BULK_CP>
        mixedLayerBulkHeatCapacity; // J K⁻¹ m⁻², from atmospheric state
    ModelArrayRef<ProtectedArray::SST> sst; // sea surface temperature, ˚C
    ModelArrayRef<ProtectedArray::TF> tf; // ocean freezing point, ˚C
    ModelArrayRef<SharedArray::DELTA_HICE> deltaHi; // New ice thickness this timestep, m

    static double minc; // Minimum sea ice concentration
    static double minh; // Minimum sea ice thickness

    void newIceFormation(size_t i, const TimestepTime&);
    void lateralIceSpread(size_t i, const TimestepTime&);
    void applyLimits(size_t i, const TimestepTime&);
    void updateWrapper(size_t i, const TimestepTime& tst)
    {
        newIceFormation(i, tst);
        lateralIceSpread(i, tst);
        applyLimits(i, tst);
    }
};

} /* namespace Nextsim */

#endif /* PHYSICS_SRC_INCLUDE_ICEGROWTH_HPP */
