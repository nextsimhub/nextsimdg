/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef HIBLERSPREAD_HPP
#define HIBLERSPREAD_HPP

#include "include/Configured.hpp"
#include "include/ILateralIceSpread.hpp"
#include "include/IceMinima.hpp"
#include "include/constants.hpp"

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

    void update(const TimestepTime& tstep) override;

private:
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
