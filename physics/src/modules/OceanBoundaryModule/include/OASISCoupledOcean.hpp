/*!
 * @file OASISCoupledOcean.hpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef OASISCOUPLEDOCEAN_HPP
#define OASISCOUPLEDOCEAN_HPP

#include "include/IOceanBoundary.hpp"
#include "include/OASISCoupled.hpp"

namespace Nextsim {

//* Ocean boundary data values that are hardcoded.
class OASISCoupledOcean : public IOceanBoundary, OASISCoupled {
public:
    OASISCoupledOcean();
    ~OASISCoupledOcean() { OASISCoupled::~OASISCoupled(); }

    std::string getName() const override { return "OASISCoupledOcean"; }
    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override;
    void setMetadata(const ModelMetadata& metadata) override;
};

} /* namespace Nextsim */

#endif /* OASISCOUPLEDOCEAN_HPP */
