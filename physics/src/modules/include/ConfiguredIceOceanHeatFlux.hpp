/*!
 * @file ConfiguredIceOceanHeatFlux.hpp
 *
 * @date 5 Jul 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGUREDICEOCEANHEATFLUX_HPP
#define CONFIGUREDICEOCEANHEATFLUX_HPP

#include "IIceOceanHeatFlux.hpp"

namespace Nextsim {

class ConfiguredIceOceanHeatFlux : public IIceOceanHeatFlux,
                                   public Configured<ConfiguredIceOceanHeatFlux> {
public:
    ConfiguredIceOceanHeatFlux();
    ~ConfiguredIceOceanHeatFlux() = default;

    enum {
        QIO_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "ConfiguredIceOceanHeatFlux"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;

    void update(const TimestepTime& tst) override { }

private:
    static double qio0;
    double qioConf;
};

} /* namespace Nextsim */

#endif /* CONFIGUREDICEOCEANHEATFLUX_HPP */
