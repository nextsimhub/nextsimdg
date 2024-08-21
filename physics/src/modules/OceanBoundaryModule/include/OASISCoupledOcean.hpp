/*!
 * @file OASISCoupledOcean.hpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef OASISCOUPLEDOCEAN_HPP
#define OASISCOUPLEDOCEAN_HPP

#include "include/IFreezingPoint.hpp"
#include "include/IOceanBoundary.hpp"
#include "include/Module.hpp"
#include "include/OASISCoupled.hpp"
#include "include/SlabOcean.hpp"

namespace Nextsim {

static const std::string moduleName = "OASISCoupledOcean";

static const std::string layerDepthConfigKey = moduleName + ".layer_depth";
static const std::string exchangeFirstLayerConfigKey = moduleName + ".exchange_first_layer";

static const double FIRST_LAYER_DEPTH = 1.; // There really is no sensible default here(!)
static const bool EXCHANGE_FIRST_LAYER = false;

static const std::string SSTKeyDefault = "I_SST";
static const std::string SSSKeyDefault = "I_SSS";
static const std::string UOceanKeyDefault = "I_Uocn";
static const std::string VOceanKeyDefault = "I_Vocn";
static const std::string SSHKeyDefault = "I_SSH";
static const std::string MLDKeyDefault = "I_MLD"; // This one is optional
static const std::string TauXKeyDefault = "I_taux";
static const std::string TauYKeyDefault = "I_tauy";
static const std::string TauModKeyDefault = "I_taumod";
static const std::string EMPKeyDefault = "I_fwflux";
static const std::string QNoSunKeyDefault = "I_rsnos";
static const std::string QSWKeyDefault = "I_rsso";
static const std::string SFluxKeyDefault = "I_sfi";
static const std::string CIceKeyDefault = "I_conc";

static const std::string SSTConfigKey = ".sea_surface_temperature";
static const std::string SSSConfigKey = ".sea_surface_salinity";
static const std::string UOceanConfigKey = ".ocean_u_velocity";
static const std::string VOceanConfigKey = ".ocean_v_velocity";
static const std::string SSHConfigKey = ".sea_surface_height";
static const std::string MLDConfigKey = ".first_ocean_layer_depth"; // This one is optional
static const std::string TauXConfigKey = ".ice_ocean_stress_x";
static const std::string TauYConfigKey = ".ice_ocean_stress_y";
static const std::string TauModConfigKey = ".ice_ocean_stress_modulo";
static const std::string EMPConfigKey = ".fresh_water_flux";
static const std::string QNoSunConfigKey = ".non_solar_heatflux";
static const std::string QSWConfigKey = ".short_wave_flux";
static const std::string SFluxConfigKey = ".salt_flux";
static const std::string CIceConfigKey = ".sea_ice_concentration";

//* Ocean boundary data values that are hardcoded.
class OASISCoupledOcean : public IOceanBoundary,
                          public OASISCoupled,
                          public Configured<OASISCoupledOcean> {
public:
    OASISCoupledOcean()
        : IOceanBoundary()
    {
    }
    ~OASISCoupledOcean() { OASISCoupled::~OASISCoupled(); }

    std::string getName() const override { return moduleName; }
    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override;
    void setMetadata(const ModelMetadata& metadata) override;

    void configure() override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

private:
    int bundleSize = 1; // Always "unbundled", as per the OASIS manual
    double firstLayerDepth = FIRST_LAYER_DEPTH;

    SlabOcean slabOcean;

    void updateTf(size_t i, const TimestepTime& tst)
    {
        tf[i] = Module::getImplementation<IFreezingPoint>()(sss[i]);
    }

    // A map to relate the strings in the namcouple file to the numbers def_var spits out
    std::map<std::string, int> couplingId;
    std::string SSTKey, SSSKey, UOceanKey, VOceanKey, SSHKey, MLDKey, TauXKey, TauYKey, TauModKey,
        EMPKey, QNoSunKey, QSWKey, SFluxKey, CIceKey;

    std::vector<std::string> cplStringsIn
        = { SSTKey, SSSKey, UOceanKey, VOceanKey, SSHKey }; // MLDKey can be added to this one
    std::vector<std::string> cplStringsOut
        = { TauXKey, TauYKey, TauModKey, EMPKey, QNoSunKey, QSWKey, SFluxKey, CIceKey };
};

} /* namespace Nextsim */

#endif /* OASISCOUPLEDOCEAN_HPP */
