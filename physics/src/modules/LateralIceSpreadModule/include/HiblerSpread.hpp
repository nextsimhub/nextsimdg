/*!
 * @file HiblerSpread.hpp
 *
 * @date 04 Jun 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef HIBLERSPREAD_HPP
#define HIBLERSPREAD_HPP

#include "include/Configured.hpp"
#include "include/ILateralIceSpread.hpp"

namespace Nextsim {

//! A class implementing the lateral spread of ice according to Hibler (1979)
class HiblerSpread : public ILateralIceSpread, public Configured<HiblerSpread> {
public:
    HiblerSpread()
        : ILateralIceSpread()
        , hice0(getStore())
        , mixedLayerBulkHeatCapacity(getStore())
        , sst(getStore())
        , tf(getStore())
    {
    }
    virtual ~HiblerSpread() = default;

    void configure() override;
    enum {
        H0_KEY,
        PHIM_KEY,
    };

    ConfigMap getConfiguration() const override;

    ModelState getStateDiagnostic() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap&, bool getAll);

    void update(const TimestepTime& tstep) override
    {
        overElements(
            [this](size_t i, const TimestepTime& tst) { this->updateWrapper(i, tst); }, tstep);
    }

private:
    void updateWrapper(size_t i, const TimestepTime& tst)
    {
        newIceFormation(i, tst);
        lateralIceSpread(i, tst);
        applyLimits(i, tst);
    }

    /*!
     * Updates the freezing of open water for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     * @param hice The ice-average ice thickness.
     * @param hsnow The ice-average snow thickness.
     * @param deltaHi The change in ice thickness this timestep.
     * @param newIce The positive change in ice thickness this timestep.
     * @param cice The ice concentration.
     * @param qow The open water heat flux.
     * @param deltaCFreeze The change in concentration due to freezing.
     */
    void freeze(const TimestepTime& tstep, double hice, double hsnow, double deltaHi, double newIce,
        double& cice, double& qow, double& deltaCfreeze);

    /*!
     * Updates the lateral melting of ice for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     * @param hice The ice-average ice thickness.
     * @param hsnow The ice-average snow thickness.
     * @param deltaHi The change in ice thickness this timestep.
     * @param cice The ice concentration.
     * @param qow The open water heat flux.
     * @param deltaCmelt The change in concentration due to melting.
     */
    void melt(const TimestepTime& tstep, double hice, double hsnow, double deltaHi, double& cice,
        double& qow, double& deltaCmelt);
    void newIceFormation(size_t i, const TimestepTime& tst);
    double updateThickness(double& thick, double newConc, double deltaC, double deltaV);
    void lateralIceSpread(size_t i, const TimestepTime& tstep);
    void applyLimits(size_t i, const TimestepTime& tstep);

    static double h0;
    static double phiM;

    ModelArrayRef<Protected::ML_BULK_CP>
        mixedLayerBulkHeatCapacity; // J K⁻¹ m⁻², from atmospheric state
    ModelArrayRef<Protected::SST> sst; // sea surface temperature, ˚C
    ModelArrayRef<Protected::TF> tf; // ocean freezing point, ˚C
    ModelArrayRef<Protected::HTRUE_ICE> hice0; // Timestep initial true ice thickness, m
};

}

#endif /* HIBLERSPREAD_HPP */
