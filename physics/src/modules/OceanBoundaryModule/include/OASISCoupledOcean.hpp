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

namespace Nextsim {

//* Ocean boundary data values that are hardcoded.
class OASISCoupledOcean : public IOceanBoundary {
public:
    OASISCoupledOcean();

    std::string getName() const override { return "OASISCoupledOcean"; }
    void setData(const ModelState::DataMap&) override;
    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override;
};

} /* namespace Nextsim */

#endif /* OASISCOUPLEDOCEAN_HPP */
