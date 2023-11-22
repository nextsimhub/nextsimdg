/*!
 * @file ConstantOceanBoundary.hpp
 *
 * @date Sep 26, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONSTANTOCEANBOUNDARY_HPP
#define CONSTANTOCEANBOUNDARY_HPP

#include "IOceanBoundary.hpp"

namespace Nextsim {

//* Ocean boundary data values that are hardcoded.
class ConstantOceanBoundary : public IOceanBoundary {
public:
    ConstantOceanBoundary();

    std::string getName() const override { return "ConstantOceanBoundary"; }
    void setData(const ModelState::DataMap&) override;
    void updateBefore(const TimestepTime& tst) override;
    void updateAfter(const TimestepTime& tst) override;
};

} /* namespace Nextsim */

#endif /* CONSTANTOCEANBOUNDARY_HPP */
