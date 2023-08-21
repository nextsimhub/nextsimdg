/*!
 * @file CCSMIceAlbedo.hpp
 *
 * @date Sep 28, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CCSMICEALBEDO_HPP
#define CCSMICEALBEDO_HPP

#include "IIceAlbedo.hpp"
#include "include/Configured.hpp"

namespace Nextsim {

//! The implementation class for the CCSM calculation of ice surface albedo.
class CCSMIceAlbedo : public IIceAlbedo, public Configured<CCSMIceAlbedo> {
public:
    /*!
     * @brief Calculates the CCSM ice surface short wave albedo.
     *
     * @param temperature The temperature of the ice surface.
     * @param snowThickness The true snow thickness on top of the ice.
     */
    std::tuple<double, double> albedo(double temperature, double snowThickness, double i0) override;

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
