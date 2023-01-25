/*!
 * @file IceGrowth.hpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICEGROWTH_HPP
#define ICEGROWTH_HPP

#include "include/Configured.hpp"
#include "include/IIceThermodynamics.hpp"
#include "include/ILateralIceSpread.hpp"
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
    };

    void configure() override;
    ConfigMap getConfiguration() const override;

    std::string getName() const override { return "IceGrowth"; }

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    std::unordered_set<std::string> hFields() const override
    {
        return { "updated_hice", "updated_cice", "updated_hsnow" };
    }
    std::unordered_set<std::string> uFields() const override { return {}; }
    std::unordered_set<std::string> vFields() const override { return {}; }
    std::unordered_set<std::string> zFields() const override { return {}; }

    void update(const TimestepTime&);

private:
    // Vertical Growth ModelComponent & Module
    std::unique_ptr<IIceThermodynamics> iVertical;
    // Lateral Growth ModuleComponent & Module
    std::unique_ptr<ILateralIceSpread> iLateral;

    // Data fields
    // Owned, shared data fields
    HField hice; // Updated true ice thickness, m
    HField cice; // Updated ice concentration
    HField hsnow; // Updated true snow thickness, m
    HField newice; // New ice over open water this timestep, m
    HField hice0; // Timestep initial true ice thickness, m
    HField hsnow0; // Timestep initial true snow thickness, m

    HField snowMelt; // Ocean to snow transfer of freshwater kg m⁻²
    // Since ILateralSpread is purely per-element, hold Δcice here
    HField deltaCIce; // Change in ice concentration
    // Owned data fields, not shared
    HField deltaCFreeze; // New ice concentration due to freezing (+ve)
    HField deltaCMelt; // Ice concentration loss due to melting (-ve)

    ModelArrayRef<ProtectedArray::H_ICE, MARConstBackingStore>
        hIceCell; // Timestep initial cell averaged ice thickness, m
    ModelArrayRef<ProtectedArray::H_SNOW, MARConstBackingStore>
        hSnowCell; // Timestep initial cell averaged snow thickness, m
    ModelArrayRef<ProtectedArray::C_ICE, MARConstBackingStore>
        cice0; // Timestep initial ice concentration
    ModelArrayRef<SharedArray::Q_OW, MARBackingStore, RW>
        qow; // open water heat flux, from FluxCalculation
    ModelArrayRef<ProtectedArray::ML_BULK_CP, MARConstBackingStore>
        mixedLayerBulkHeatCapacity; // J K⁻¹ m⁻², from atmospheric state
    ModelArrayRef<ProtectedArray::SST, MARConstBackingStore> sst; // sea surface temperature, ˚C
    ModelArrayRef<ProtectedArray::TF, MARConstBackingStore> tf; // ocean freezing point, ˚C
    ModelArrayRef<SharedArray::DELTA_HICE, MARBackingStore>
        deltaHi; // New ice thickness this timestep, m

    void newIceFormation(size_t i, const TimestepTime&);
    void lateralIceSpread(size_t i, const TimestepTime&);
    void applyLimits(size_t i, const TimestepTime&);
    void updateWrapper(size_t i, const TimestepTime& tst)
    {
        newIceFormation(i, tst);
        lateralIceSpread(i, tst);
        applyLimits(i, tst);
    }
    void initializeThicknesses(size_t i, const TimestepTime&);
};

} /* namespace Nextsim */

#endif /* ICEGROWTH_HPP */
