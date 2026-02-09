/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef SRC_INCLUDE_SMU2ICEALBEDO_HPP
#define SRC_INCLUDE_SMU2ICEALBEDO_HPP

#include "include/IIceAlbedo.hpp"

namespace Nextsim {

//! The implementation class for the SMU calculation of ice surface albedo
// with variable snow albedo.
class SMU2IceAlbedo : public IIceAlbedo, public Configured<SMU2IceAlbedo> {
public:
    std::string getName() const override;

    void configure() override;
    ConfigMap getConfiguration() const override;
    void update(const TimestepTime& tst) override;

private:
    static double i0;
};

}

#endif /* SRC_INCLUDE_SMU2ICEALBEDO_HPP */
