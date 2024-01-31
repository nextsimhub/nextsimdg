/*!
 * @file CGDynamicsKernel.hpp
 *
 * @date Jan 31, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CGDYNAMICSKERNEL_HPP
#define CGDYNAMICSKERNEL_HPP

#include "DynamicsKernel.hpp"

// Import this from the build system *somehow*
static const int CGdegree = 2;
static const int DGstressDegree = CG2DGSTRESS(CGdegree);
static const int nGauss = CGdegree + 1;
static const int CGdof = nGauss * nGauss;

namespace Nextsim {

template <int DGadvection>
class CGDynamicsKernel: public DynamicsKernel<DGadvection, 1> {
    CGDynamicsKernel() = default;
    virtual ~CGDynamicsKernel() = default;
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void setData(const std::string& name, const ModelArray& data) override;
    ModelArray getDG0Data(const std::string& name) override;
    void update(const TimestepTime& tst) override;
    void updateMomentum(const TimestepTime& tst) override;
    void prepareIteration(const DataMap& data) override;
    void projectVelocityToStrain() override;
    void calculateStressDivergence(const double scale) override;
    void applyBoundaries() override;
    void prepareAdvection() override;
    void setVelocityData(const std::string& velocityFieldName, const ModelArray& data) override;
    ModelArray getVelocityDG0Data(const std::string& name) override;
};

} /* namespace Nextsim */

#endif /* CGDYNAMICSKERNEL_HPP */
