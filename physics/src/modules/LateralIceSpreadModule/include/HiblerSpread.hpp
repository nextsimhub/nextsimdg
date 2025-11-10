/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef HIBLERSPREAD_HPP
#define HIBLERSPREAD_HPP

#include "include/Configured.hpp"
#include "include/ILateralIceSpread.hpp"
#include "include/constants.hpp"
#include "include/IceMinima.hpp"

namespace Nextsim {

//! A class implementing the lateral spread of ice according to Hibler (1979)
class HiblerSpread : public ILateralIceSpread, public Configured<HiblerSpread> {
public:
    HiblerSpread()
        : ILateralIceSpread()
        , hiceAccessor(getStore())
        , mixedLayerBulkHeatCapacityAccessor(getStore())
        , sstAccessor(getStore())
        , tfAccessor(getStore())
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
        AdvectedField& hsnow = hsnowAccessor.getHostRW();
        HField& qow = qowAccessor.getHostRW();
        AdvectedField& cice = ciceAccessor.getHostRW();
        HField& newice = newiceAccessor.getHostRW();
        AdvectedField& hice = hiceAccessor.getHostRW();
        HField& deltaCIce = deltaCIceAccessor.getHostRW();
        const HField& mixedLayerBulkHeatCapacity = mixedLayerBulkHeatCapacityAccessor.getHostRO();
        const HField& deltaHi = deltaHiAccessor.getHostRO();
        const HField& tf = tfAccessor.getHostRO();
        const HField& sst = sstAccessor.getHostRO();

        overElements(
            [&](size_t i, const TimestepTime& tst) {
                // newIceFormation
                // Flux cooling the ocean from open water
                // TODO Add assimilation fluxes here
                double coolingFlux = qow[i];
                // Temperature change of the mixed layer during this timestep
                double deltaTml = -coolingFlux / mixedLayerBulkHeatCapacity[i] * tst.step;
                // Initial temperature
                double t0 = sst[i];
                // Freezing point temperature
                double tf0 = tf[i];
                // Final temperature
                double t1 = t0 + deltaTml;

                // deal with cooling below the freezing point
                if (t1 < tf0) {
                    // Heat lost cooling the mixed layer to freezing point
                    double sensibleFlux = (tf0 - t0) / deltaTml * coolingFlux;
                    // Any heat beyond that is latent heat forming new ice
                    double latentFlux = coolingFlux - sensibleFlux;

                    qow[i] = sensibleFlux;
                    newice[i] = latentFlux * tst.step * (1 - cice[i]) / (Ice::Lf * Ice::rho);
                } else {
                    newice[i] = 0;
                }

                // lateralIceSpread
                const double deltaCMelt = melt(deltaHi[i], cice[i], hice[i]);
                const double deltaCFreeze = freeze(newice[i]);

                deltaCIce[i] = deltaCFreeze + deltaCMelt;
                cice[i] = (hice[i] > 0 || newice[i] > 0) ? cice[i] + deltaCIce[i] : 0;
                if (cice[i] >= IceMinima::c()) {
                    // The updated ice thickness must conserve volume
                    hice[i] += newice[i];
                    if (deltaCIce[i] < 0) {
                        /* Snow is lost if the concentration decreases, and energy is returned
                         * to the ocean. We reduce the snow volume by a "slice" of snow with the
                         * dimensions hs * deltaCIce. */
                        const double hs = hsnow[i] / (cice[i] - deltaCIce[i]);
                        qow[i] -= deltaCIce[i] * hs * Water::Lf * Ice::rhoSnow / tstep.step;
                        hsnow[i] += hs * deltaCIce[i];
                    } // else: Snow volume is conserved, so no change to hsnow[i]
                }

                // applyLimits
                if (cice[i] < IceMinima::c() || hice[i] < IceMinima::h()) {
                    qow[i]
                        += Water::Lf * (hice[i] * Ice::rho + hsnow[i] * Ice::rhoSnow) / tstep.step;
                    hice[i] = 0;
                    cice[i] = 0;
                    hsnow[i] = 0;
                }
            },
            tstep);
    }

private:
    /*!
     * Updates the freezing of open water for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     * @param newIce The positive change in ice thickness this timestep.
     * @param deltaCFreeze The change in concentration due to freezing.
     */
    double freeze(double newIce);

    /*!
     * Updates the lateral melting of ice for the timestep.
     *
     * @param tStep The object containing the timestep start and duration times.
     * @param deltaHi The change in ice thickness this timestep.
     * @param cice The ice concentration.
     * @param hice The ice-average ice thickness.
     */
    double melt(double deltaHi, double cice, double hice);
    void newIceFormation(size_t i, const TimestepTime& tst);
    void lateralIceSpread(size_t i, const TimestepTime& tstep);
    void applyLimits(size_t i, const TimestepTime& tstep);

    static double h0;
    static double phiM;

    ModelArrayAccessor<Protected::ML_BULK_CP>
        mixedLayerBulkHeatCapacityAccessor; // J K⁻¹ m⁻², from atmospheric state
    ModelArrayAccessor<Protected::SST> sstAccessor; // sea surface temperature, ˚C
    ModelArrayAccessor<Protected::TF> tfAccessor; // ocean freezing point, ˚C
    ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor; // Timestep initial true ice thickness, m
};

}

#endif /* HIBLERSPREAD_HPP */
