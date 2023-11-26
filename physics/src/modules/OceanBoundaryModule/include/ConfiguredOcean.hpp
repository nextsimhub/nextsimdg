/*!
 * @file ConfiguredOcean.hpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
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

    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override;

private:
    static double sst0;
    static double sss0;
    static double mld0;
    static double u0;
    static double v0;

    // External SS* fields to feed the slab ocean
    HField sstExt;
    HField sssExt;

    // We need a slab ocean in this implementation
    SlabOcean slabOcean;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDOCEAN_HPP */
