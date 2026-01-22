/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_SMUICEALBEDO_HPP
#define SRC_INCLUDE_SMUICEALBEDO_HPP

#include "include/Configured.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

//! The implementation class for the SMU calculation of ice surface albedo
// with constant snow albedo.
class SMUIceAlbedo : public IIceAlbedo, public Configured<SMUIceAlbedo> {
public:
    void configure() override;
    ConfigMap getConfiguration() const override;
    void update(const TimestepTime& tst) override;

    enum {
        I0_KEY,
    };

    /*!
     * @brief Calculates the SMU ice surface short wave albedo with constant
     * snow albedo.
     *
     * @param temperature The temperature of the ice surface.
     * @param snowThickness The true snow thickness on top of the ice.
     * @param i0 The fraction of short-wave radiation that can penetrate bare ice (not taking snow
     * cover into account).
     */
    std::tuple<double, double> surfaceShortWaveBalance(
        double temperature, double snowThickness, double i0);

private:
    static double m_i0;
};

}

#endif /* SRC_INCLUDE_SMUICEALBEDO_HPP */
