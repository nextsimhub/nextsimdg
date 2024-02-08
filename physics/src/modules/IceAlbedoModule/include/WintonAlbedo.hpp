/*!
 * @file WintonAlbedo.hpp
 *
 * @date Sep 11, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef WINTONALBEDO_HPP
#define WINTONALBEDO_HPP

#include "include/Configured.hpp"
#include "include/IIceAlbedo.hpp"

namespace Nextsim {

//! The implementation class for the Winton (2000) calculation of ice surface albedo.
class WintonAlbedo : public IIceAlbedo, public Configured<WintonAlbedo> {
public:
    /*!
     * @brief Calculates the ice surface short wave albedo and fraction of penetrating
     * short-wave radiation following Winton (2000).
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
    static double meltAlbedo;
};

}

#endif /* WINTONALBEDO_HPP */
