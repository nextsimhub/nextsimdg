/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef CCSMICEALBEDO_HPP
#define CCSMICEALBEDO_HPP

#include "include/Configured.hpp"
#include "include/IIceAlbedo.hpp"

namespace Nextsim {

//! The implementation class for the CCSM calculation of ice surface albedo.
class CCSMIceAlbedo : public IIceAlbedo, public Configured<CCSMIceAlbedo> {
public:
    std::string getName() const override;

    void configure() override;
    ConfigMap getConfiguration() const override;

    void update(const TimestepTime& tst) override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

private:
    static double iceAlbedo;
    static double snowAlbedo;
    static double i0;
};

}

#endif /* CCSMICEALBEDO_HPP */
