/*!
 * @file ConfiguredOcean.hpp
 *
 * @date Aug 31, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDOCEAN_HPP
#define CONFIGUREDOCEAN_HPP

#include "OceanState.hpp"

#include "include/Configured.hpp"

namespace Nextsim {

class ConfiguredOcean : public OceanState, public Configured<ConfiguredOcean> {
public:
    ConfiguredOcean() = default;
    ~ConfiguredOcean() = default;

    enum {
        SST_KEY,
        SSS_KEY,
        MLD_KEY,
    };

    void setData(const ModelState::DataMap&) override { }
    std::string getName() const override { return "ConfiguredOcean"; }

    void configure() override;

protected:
    //! Performs the implementation specific updates. Does nothing.
    void updateSpecial(const TimestepTime&) override { }

private:
    static double sst0;
    static double sss0;
    static double mld0;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDOCEAN_HPP */
