/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FLUXCONFIGUREDOCEAN_HPP
#define FLUXCONFIGUREDOCEAN_HPP

#include "include/IOceanBoundary.hpp"

#include "include/Configured.hpp"

namespace Nextsim {

//! A class to provide constant oceanic forcings that can be configured at run time.
class FluxConfiguredOcean : public IOceanBoundary, public Configured<FluxConfiguredOcean> {
public:
    FluxConfiguredOcean() = default;
    ~FluxConfiguredOcean() = default;

    enum {
        QIO_KEY,
        SST_KEY,
        SSS_KEY,
        MLD_KEY,
        CURRENTU_KEY,
        CURRENTV_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "FluxConfiguredOcean"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;
    ConfigMap getConfiguration() const override;

    void updateBefore(const TimestepTime& tst) override { }
    void updateAfter(const TimestepTime& tst) override { }

private:
    static FloatType qio0;
    static FloatType sst0;
    static FloatType sss0;
    static FloatType mld0;
    static FloatType u0;
    static FloatType v0;
};

} /* namespace Nextsim */

#endif /* FLUXCONFIGUREDOCEAN_HPP */
