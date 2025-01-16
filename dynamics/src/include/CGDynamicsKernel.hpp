/*!
 * @file CGDynamicsKernel.hpp
 *
 * @date 06 Dec 2024
 * @author Tim Spain <timothy.spain@nersc.no>
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

public:
    CGDynamicsKernel() { }
    virtual ~CGDynamicsKernel() = default;
    void initialise(const ModelArray& coords, bool isSpherical, const ModelArray& mask) override;

    void setData(const std::string& name, const ModelArray& data) override;
    ModelArray getDG0Data(const std::string& name) const override;
    void ComputeGradientOfSeaSurfaceHeight(const DGVector<1>& seaSurfaceHeight);
    void prepareIteration(const DataMap& data) override;
    void projectVelocityToStrain() override;
    void stressDivergence() override;
    void applyBoundaries() override;
    void prepareAdvection() override;

    virtual inline double getIceOceanStressElement(const std::string& name, const int i) const = 0;
    CGVector<CGdegree> getIceOceanStress(const std::string& name) const
    {
        if (name != uIOStressName && name != vIOStressName)
            throw std::logic_error(std::string(__func__) + " called with an unknown argument "
                + name + ". Only " + uIOStressName + " and " + vIOStressName + " are supported\n");

        CGVector<CGdegree> tau;
        tau.resizeLike(u);

#pragma omp parallel for
        for (int i = 0; i < tau.rows(); ++i)
            tau(i) = getIceOceanStressElement(name, i);

        return tau;
    }

protected:
    void addStressTensorCell(const size_t eid, const size_t cx, const size_t cy);
    void dirichletZero(CGVector<CGdegree>&) const;
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

    std::unique_ptr<ParametricMomentumMap<CGdegree, DGadvection>> pmap;
};

} /* namespace Nextsim */

#endif /* CGDYNAMICSKERNEL_HPP */
