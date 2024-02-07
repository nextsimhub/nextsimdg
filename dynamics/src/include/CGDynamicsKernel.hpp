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
class CGDynamicsKernel : public DynamicsKernel<DGadvection, DGstressDegree> {
    using DynamicsKernel<DGadvection, DGstressDegree>::smesh;
public:
    CGDynamicsKernel() = default;
    virtual ~CGDynamicsKernel() = default;
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void setData(const std::string& name, const ModelArray& data) override;
    ModelArray getDG0Data(const std::string& name) override;
    void update(const TimestepTime& tst) override;
    void updateMomentum(const TimestepTime& tst) override;
    void prepareIteration(const typename DynamicsKernel<DGadvection, DGstressDegree>::DataMap& data) override;
    void projectVelocityToStrain() override;
    void stressDivergence(const double scale) override;
    void applyBoundaries() override;
    void prepareAdvection() override;
protected:
    void addStressTensorCell(const double scale, const size_t eid, const size_t cx, const size_t cy);
    void dirichletZero(CGVector<CGdegree>&) const;
    // CG ice velocity
    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    // CG ice thickness and concentration
    CGVector<CGdegree> cgA;
    CGVector<CGdegree> cgH;

    // divergence of stress
    CGVector<CGdegree> dStressX;
    CGVector<CGdegree> dStressY;

    // Ocean velocity
    CGVector<CGdegree> uOcean;
    CGVector<CGdegree> vOcean;

    // Atmospheric wind velocity
    CGVector<CGdegree> uAtmos;
    CGVector<CGdegree> vAtmos;

//    CGParametricMomentum<CGdegree>* momentum;

};

} /* namespace Nextsim */

#endif /* CGDYNAMICSKERNEL_HPP */
