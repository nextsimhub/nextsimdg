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
public:
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
protected:
    void addStressTensorCell(const double scale, const size_t eid, const size_t cx, const size_t cy);

    // CG ice velocity
    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    // CG ice thickness and concentration
    CGVector<CGdegree> cgA;
    CGVector<CGdegree> cgH;

    // divergence of stress
    CGVector<CGdegree> dStressX;
    CGVector<CGdegree> dStressY;

    CGParametricMomentum<CGdegree>* momentum;

    ParametricMomentumMap<CGdegree> pmap;
};

} /* namespace Nextsim */

#endif /* CGDYNAMICSKERNEL_HPP */
