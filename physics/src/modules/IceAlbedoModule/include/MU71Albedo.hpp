/*!
 * @author  Einar Örn Ólason <einar.olason@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef SEASONALICEALBEDO_HPP
#define SEASONALICEALBEDO_HPP

#include "include/IIceAlbedo.hpp"
#include "include/MonthlyCubicBSpline.hpp"

namespace Nextsim {

/*!
 * @brief The implementation class for the a seasonal albedo following Maykut and Untersteiener's
 * (1971) table 1. Only useful for comparison with that paper and derived setups.
 */
class MU71Albedo : public IIceAlbedo, public Configured<MU71Albedo> {

public:
    MU71Albedo();

    std::string getName() const override;
    void configure() override;
    ConfigMap getConfiguration() const override;
    void update(const TimestepTime& tst) override;

private:
    monthlyCubicBSpline snowAlbedo;

    static FloatType i0;
};

}

#endif /* SEASONALICEALBEDO_HPP */
