/*!
 * @file CGDynamicsKernel.hpp
 *
 * @date 27 Aug 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef CGDYNAMICSKERNEL_HPP
#define CGDYNAMICSKERNEL_HPP

#include "DynamicsKernel.hpp"

#ifndef CGDEGREE
#define CGDEGREE 2
#define DGSTRESSCOMP (CG2DGSTRESS(CGDEGREE))
#endif

static const int CGdegree = CGDEGREE;
static const int DGstressComp = DGSTRESSCOMP;
static const int nGauss = CGdegree + 1;
static const int CGdof = nGauss * nGauss;

namespace Nextsim {

template <int DGadvection>
class CGDynamicsKernel : public DynamicsKernel<DGadvection, DGstressComp> {
protected:
    using DynamicsKernel<DGadvection, DGstressComp>::s11;
    using DynamicsKernel<DGadvection, DGstressComp>::s12;
    using DynamicsKernel<DGadvection, DGstressComp>::s22;
    using DynamicsKernel<DGadvection, DGstressComp>::e11;
    using DynamicsKernel<DGadvection, DGstressComp>::e12;
    using DynamicsKernel<DGadvection, DGstressComp>::e22;
    using DynamicsKernel<DGadvection, DGstressComp>::smesh;
    using DynamicsKernel<DGadvection, DGstressComp>::dgtransport;
    using typename DynamicsKernel<DGadvection, DGstressComp>::DataMap;
    using DynamicsKernel<DGadvection, DGstressComp>::cosOceanAngle;
    using DynamicsKernel<DGadvection, DGstressComp>::sinOceanAngle;

public:
    CGDynamicsKernel(double cosOceanAngleIn, double sinOceanAngleIn,
        const DynamicsParameters& paramsIn, CGVector<CGdegree>& uStressRef,
        CGVector<CGdegree>& vStressRef)
        : DynamicsKernel<DGadvection, DGstressComp>(cosOceanAngleIn, sinOceanAngleIn, paramsIn)
        , pmap(nullptr)
        , uStress(uStressRef)
        , vStress(vStressRef)
    {
    }
    CGDynamicsKernel(
        double cosOceanAngleIn, const DynamicsParameters& paramsIn, double sinOceanAngleIn)
        : CGDynamicsKernel(cosOceanAngleIn, sinOceanAngleIn, paramsIn, u, v)
    {
    }
    virtual ~CGDynamicsKernel() = default;
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;
    void setData(const std::string& name, const ModelArray& data) override;
    ModelArray getDG0Data(const std::string& name) const override;
    void prepareIteration(const DataMap& data) override;
    void projectVelocityToStrain() override;
    void stressDivergence() override;
    void applyBoundaries() override;
    void prepareAdvection() override;

    virtual std::pair<CGVector<CGdegree>, CGVector<CGdegree>> getIceOceanStress() const;

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

    // Velocities used to calculate ice-ocean stress. The references must be set in derived classes.
    CGVector<CGdegree>& uStress;
    CGVector<CGdegree>& vStress;

    ParametricMomentumMap<CGdegree>* pmap;
};

} /* namespace Nextsim */

#endif /* CGDYNAMICSKERNEL_HPP */
