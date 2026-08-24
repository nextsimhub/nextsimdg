/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDOCEAN_HPP
#define CONFIGUREDOCEAN_HPP

#include "include/IOceanBoundary.hpp"

#include "include/Configured.hpp"
#include "include/SlabOcean.hpp"

namespace Nextsim {

//! A class to provide constant oceanic forcings that can be configured at run
//! time as physical variables.
class ConfiguredOcean : public IOceanBoundary, public Configured<ConfiguredOcean> {
public:
    ConfiguredOcean();
    ~ConfiguredOcean() = default;

    enum {
        SST_KEY,
        SSS_KEY,
        MLD_KEY,
        CURRENTU_KEY,
        CURRENTV_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "ConfiguredOcean"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;
    ConfigMap getConfiguration() const override;

    ModelState getStatePrognostic() const override;
    ModelState getStateDiagnostic() const override;
    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override;

private:
    static FloatType sst0;
    static FloatType sss0;
    static FloatType mld0;
    static FloatType u0;
    static FloatType v0;

    // External SS* fields to feed the slab ocean
    ModelArrayAccessor<Protected::EXT_SST, RW> sstExtAccessor;
    ModelArrayAccessor<Protected::EXT_SSS, RW> sssExtAccessor;

    ModelArrayAccessor<Protected::SLAB_SST, RO> sstSlabAccessor;
    ModelArrayAccessor<Protected::SLAB_SSS, RO> sssSlabAccessor;

    // We need a slab ocean in this implementation
    SlabOcean slabOcean;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDOCEAN_HPP */
