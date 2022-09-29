/*!
 * @file ConfiguredOcean.hpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDOCEAN_HPP
#define CONFIGUREDOCEAN_HPP

#include "IOceanBoundary.hpp"

#include "include/Configured.hpp"

namespace Nextsim {

//! A class to provide constant oceanic forcings that can be configured at run time.
class ConfiguredOcean : public IOceanBoundary, public Configured<ConfiguredOcean> {
public:
    ConfiguredOcean() = default;
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

    void updateBefore(const TimestepTime& tst) override { }
    void updateAfter(const TimestepTime& tst) override { }

private:
    static double qio;
    static double sst0;
    static double sss0;
    static double mld0;
    static double u0;
    static double v0;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDOCEAN_HPP */
