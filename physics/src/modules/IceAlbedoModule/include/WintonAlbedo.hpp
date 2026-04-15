/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Einar Ólason <einar.olason@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef WINTONALBEDO_HPP
#define WINTONALBEDO_HPP

#include "include/Configured.hpp"
#include "include/IIceAlbedo.hpp"

namespace Nextsim {

//! The implementation class for the Winton (2000) calculation of ice surface albedo.
class WintonAlbedo : public IIceAlbedo, public Configured<WintonAlbedo> {
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
    static double meltAlbedo;
    static double i0;
};

}

#endif /* WINTONALBEDO_HPP */
