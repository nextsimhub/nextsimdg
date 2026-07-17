/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
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
    std::string getName() const override;
    void configure() override;
    ConfigMap getConfiguration() const override;
    void update(const TimestepTime& tst) override;

private:
    static FloatType i0;
};

}

#endif /* SRC_INCLUDE_SMUICEALBEDO_HPP */
