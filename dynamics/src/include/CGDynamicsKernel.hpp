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
static constexpr int CGdegree = 1;
static constexpr int DGstressDegree = CG2DGSTRESS(CGdegree);
static constexpr int nGauss = CGdegree + 1;
static constexpr int CGdof = nGauss * nGauss;

// for testing only
// define to use the Kokkos MEVP (CPU or GPU depending on cmake options) kernel, otherwise the OpenMP CPU kernel is used
#define USE_KOKKOS_KERNEL
// use Kokkos 1D or 2D parallel loop for strain and divergence computation
#define LOOP_1D
// number of dG degrees of freedom for advection; to change the cG degree modify CGdegree above
static constexpr int DGdof = 3;
// measure individual operations; should be disabled to measure the full MEVP iteration because the timers introduce 
// additional synchronization that causes a measurable slowdown
constexpr bool MEASURE_DETAILED = true;

namespace Nextsim {

template <int DGadvection>
class CGDynamicsKernel : public DynamicsKernel<DGadvection, DGstressDegree> {
protected:
    using DynamicsKernel<DGadvection, DGstressDegree>::s11;
    using DynamicsKernel<DGadvection, DGstressDegree>::s12;
    using DynamicsKernel<DGadvection, DGstressDegree>::s22;
    using DynamicsKernel<DGadvection, DGstressDegree>::e11;
    using DynamicsKernel<DGadvection, DGstressDegree>::e12;
    using DynamicsKernel<DGadvection, DGstressDegree>::e22;
    using DynamicsKernel<DGadvection, DGstressDegree>::smesh;
    using DynamicsKernel<DGadvection, DGstressDegree>::dgtransport;
    using typename DynamicsKernel<DGadvection, DGstressDegree>::DataMap;

public:
    CGDynamicsKernel()
        : pmap(nullptr)
    {
    }
    virtual ~CGDynamicsKernel() = default;
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void setData(const std::string& name, const ModelArray& data) override;
    ModelArray getDG0Data(const std::string& name) override;
    void prepareIteration(const DataMap& data) override;
    void projectVelocityToStrain() override;
    void stressDivergence() override;
    void applyBoundaries() override;
    void prepareAdvection() override;

protected:
    void addStressTensorCell(const size_t eid, const size_t cx, const size_t cy);
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

    ParametricMomentumMap<CGdegree>* pmap;
};

} /* namespace Nextsim */

#endif /* CGDYNAMICSKERNEL_HPP */
