/*!
 * @file CCSMIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CCSMICEALBEDO_HPP
#define CCSMICEALBEDO_HPP

#include "include/Configured.hpp"
#include "include/IIceAlbedo.hpp"

namespace Nextsim {

//! The implementation class for the CCSM calculation of ice surface albedo.
class CCSMIceAlbedo : public IIceAlbedo, public Configured<CCSMIceAlbedo> {
public:
    /*!
     * @brief Calculates the CCSM3 ice surface short wave albedo and fraction of penetrating
     * short-wave radiation.
     * @param temperature The temperature of the ice surface.
     * @param snowThickness The true snow thickness on top of the ice.
     * @param i0 The fraction of short-wave radiation that can penetrate bare ice (not taking snow
     * cover into account).
     */
    std::tuple<double, double> surfaceShortWaveBalance(
        double temperature, double snowThickness, double i0) override;

    void configure() override;

    ConfigMap getConfiguration() const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

private:
    static double iceAlbedo;
    static double snowAlbedo;
};

}

#endif /* CCSMICEALBEDO_HPP */
