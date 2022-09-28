/*!
 * @file ConstantAtmosphereBoundary.hpp
 *
 * @date Sep 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONSTANTATMOSPHEREBOUNDARY_HPP
#define CONSTANTATMOSPHEREBOUNDARY_HPP

#include "include/IAtmosphereBoundary.hpp"

namespace Nextsim {

//* A set of atmospheric inputs that are fixed at compile time.
class ConstantAtmosphereBoundary : public IAtmosphereBoundary {
public:
    ConstantAtmosphereBoundary();

    std::string getName() const override { return "ConstantAtmosphereBoundary"; }
    void setData(const ModelState::DataMap&) override;
    void update(const TimestepTime& tst) override;
};

} /* namespace Nextsim */

#endif /* CONSTANTATMOSPHEREBOUNDARY_HPP */
