/*!
 * @file CGDynamicsKernel.hpp
 *
 * @date 27 Mar 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef CGDYNAMICSKERNEL_HPP
#define CGDYNAMICSKERNEL_HPP

#include "DynamicsKernel.hpp"
#include "DynamicsParameters.hpp"

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

public:
    CGDynamicsKernel(const DynamicsParameters& params);
    virtual ~CGDynamicsKernel() = default;
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    void setData(const std::string& name, const ModelArray& data) override;
    ModelArray getDG0Data(const std::string& name) const override;
    void computeGradientOfSeaSurfaceHeight(const DGVector<1>& seaSurfaceHeight);
    void prepareIteration(const DataMap& data) override;
    void projectVelocityToStrain() override;
    void stressDivergence() override;
    void applyBoundaries() override;
    void prepareAdvection() override;

protected:
    void addStressTensorCell(const size_t eid, const size_t cx, const size_t cy);
    void dirichletZero(CGVector<CGdegree>&) const;
    void updateIceOceanStress(const CGVector<CGdegree>& uIce, const CGVector<CGdegree>& vIce);

    // CG ice velocity
    CGVector<CGdegree> u;
    CGVector<CGdegree> v;

    // CG ice thickness and concentration
    CGVector<CGdegree> cgA;
    CGVector<CGdegree> cgH;

    // CG gradient of the seaSurfaceHeight
    CGVector<CGdegree> xGradSeaSurfaceHeight;
    CGVector<CGdegree> yGradSeaSurfaceHeight;

    // divergence of stress
    CGVector<CGdegree> dStressX;
    CGVector<CGdegree> dStressY;

    // Ocean velocity
    CGVector<CGdegree> uOcean;
    CGVector<CGdegree> vOcean;

    // Atmospheric wind velocity
    CGVector<CGdegree> uAtmos;
    CGVector<CGdegree> vAtmos;

    // ice ocean stresses
    CGVector<CGdegree> uIceOceanStress;
    CGVector<CGdegree> vIceOceanStress;

    double cosOceanAngle = std::numeric_limits<double>::quiet_NaN();
    double sinOceanAngle = std::numeric_limits<double>::quiet_NaN();
    const DynamicsParameters& baseParams;

    std::unique_ptr<ParametricMomentumMap<CGdegree, DGadvection>> pmap;

    CGVector<CGdegree>& ma2cg(const ModelArray& maData, CGVector<CGdegree>& cgData)
    {
        DGVector<DGadvection> dgtmp(*smesh);
        dgtmp.zero();
        DGModelArray::ma2dg(maData, dgtmp);
        Nextsim::Interpolations::DG2CG(*smesh, cgData, dgtmp);
        return cgData;
    }
};

} /* namespace Nextsim */

#endif /* CGDYNAMICSKERNEL_HPP */
