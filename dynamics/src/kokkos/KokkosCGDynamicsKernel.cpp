/*!
 * @file KokkosCGDynamicsKernel.cpp
 * @date 02 Jun 2025
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosCGDynamicsKernel.hpp"
#include "include/KokkosDGLimit.hpp"
#include "include/KokkosTimer.hpp"

namespace Nextsim {
/*************************************************************/
template <int DGadvection>
KokkosCGDynamicsKernel<DGadvection>::KokkosCGDynamicsKernel(const DynamicsParameters& params)
    : CGDynamicsKernel<DGadvection>(params)
{
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::initialise(
    const ModelArray& coords, bool isSpherical, const ModelArray& mask)
{
    CGDynamicsKernel<DGadvection>::initialise(coords, isSpherical, mask);

    // velocity
    std::tie(uHost, uDevice) = makeKokkosDualView("u", this->u);
    std::tie(vHost, vDevice) = makeKokkosDualView("v", this->v);

    std::tie(cgHHost, cgHDevice) = makeKokkosDualView("cgH", this->cgH);
    std::tie(cgAHost, cgADevice) = makeKokkosDualView("cgA", this->cgA);

    std::tie(xGradSeaSurfaceHeightHost, xGradSeaSurfaceHeightDevice)
        = makeKokkosDualView("xGradSeaSurfaceHeight", this->xGradSeaSurfaceHeight);
    std::tie(yGradSeaSurfaceHeightHost, yGradSeaSurfaceHeightDevice)
        = makeKokkosDualView("yGradSeaSurfaceHeight", this->yGradSeaSurfaceHeight);

    std::tie(dStressXHost, dStressXDevice) = makeKokkosDualView("dStressX", this->dStressX);
    std::tie(dStressYHost, dStressYDevice) = makeKokkosDualView("dStressY", this->dStressY);

    std::tie(uOceanHost, uOceanDevice) = makeKokkosDualView("uOcean", this->uOcean);
    std::tie(vOceanHost, vOceanDevice) = makeKokkosDualView("vOcean", this->vOcean);

    std::tie(uAtmosHost, uAtmosDevice) = makeKokkosDualView("uAtmos", this->uAtmos);
    std::tie(vAtmosHost, vAtmosDevice) = makeKokkosDualView("vAtmos", this->vAtmos);

    std::tie(uIceOceanStressHost, uIceOceanStressDevice)
        = makeKokkosDualView("uIceOceanStress", this->uIceOceanStress);
    std::tie(vIceOceanStressHost, vIceOceanStressDevice)
        = makeKokkosDualView("vIceOceanStress", this->vIceOceanStress);

    // stress
    std::tie(s11Host, s11Device) = makeKokkosDualView("s11", this->s11);
    std::tie(s12Host, s12Device) = makeKokkosDualView("s12", this->s12);
    std::tie(s22Host, s22Device) = makeKokkosDualView("s22", this->s22);
    std::tie(e11Host, e11Device) = makeKokkosDualView("e11", this->e11);
    std::tie(e12Host, e12Device) = makeKokkosDualView("e12", this->e12);
    std::tie(e22Host, e22Device) = makeKokkosDualView("e22", this->e22);

    tempDataAdvect.resize_by_mesh(*this->smesh);
    std::tie(tempDataHost, tempDataDevice) = makeKokkosDualView("tempData", (this->tempData));
    tempData.resize_by_mesh(*this->smesh);
    std::tie(tempDataAdvectHost, tempDataAdvectDevice)
        = makeKokkosDualView("tempDataAdvect", (this->tempDataAdvect));

    assert(this->pmap);
    divS1Device = makeKokkosDeviceViewMap("divS1", this->pmap->divS1, MakeViewOptions::DEVICE_COPY);
    divS2Device = makeKokkosDeviceViewMap("divS2", this->pmap->divS2, MakeViewOptions::DEVICE_COPY);
    divMDevice = makeKokkosDeviceViewMap("divM", this->pmap->divM, MakeViewOptions::DEVICE_COPY);
    iMgradXDevice
        = makeKokkosDeviceViewMap("iMgradX", this->pmap->iMgradX, MakeViewOptions::DEVICE_COPY);
    iMgradYDevice
        = makeKokkosDeviceViewMap("iMgradY", this->pmap->iMgradY, MakeViewOptions::DEVICE_COPY);
    iMMDevice = makeKokkosDeviceViewMap("iMM", this->pmap->iMM, MakeViewOptions::DEVICE_COPY);

    // needed for stress and momentum
    std::tie(hiceHost, hiceDevice)
        = makeKokkosDualView("hice", static_cast<DGVector<DGadvection>&>(this->hice));
    std::tie(ciceHost, ciceDevice)
        = makeKokkosDualView("cice", static_cast<DGVector<DGadvection>&>(this->cice));

    PSIAdvectDevice = makeKokkosDeviceView(
        "PSI<DGadvection, NGP>", PSI<DGadvection, NGP>, MakeViewOptions::DEVICE_COPY);
    PSIStressDevice = makeKokkosDeviceView(
        "PSI<DGstress, NGP>", PSI<DGstressComp, NGP>, MakeViewOptions::DEVICE_COPY);

    lumpedCGMassDevice = makeKokkosDeviceView(
        "lumpedcgmass", this->pmap->lumpedcgmass, MakeViewOptions::DEVICE_COPY);
    iMJwPSIDevice
        = makeKokkosDeviceViewMap("iMJwPSI", this->pmap->iMJwPSI, MakeViewOptions::DEVICE_COPY);
    iMJwPSIAdvectDevice = makeKokkosDeviceViewMap(
        "iMJwPSIAdvect", this->pmap->iMJwPSI_dam, MakeViewOptions::DEVICE_COPY);

    const size_t n_cg = static_cast<size_t>(this->pmap->cglandmask.rows());
    std::vector<bool> cgLandMask(n_cg, false);
    for (size_t i = 0; i < n_cg; ++i) {
        const auto v = this->pmap->cglandmask(i);
        // landmask is not stored as bool
        assert(v == 0 || v == 1);
        cgLandMask[i] = v != 0.0;
    }
    cgLandMaskDevice = makeKokkosDeviceBitset(cgLandMask);

    assert(this->smesh);
    meshData = std::make_unique<KokkosMesh>(*this->smesh);
    cG2DGAdvectInterpolator
        = std::make_unique<Interpolations::KokkosCG2DGInterpolator<DGadvection, CGdegree>>(
            *this->smesh);
    dG2CGAdvectInterpolator
        = std::make_unique<Interpolations::KokkosDG2CGInterpolator<CGdegree, DGadvection>>(
            *this->smesh);
    dGTransportDevice = std::make_unique<KokkosDGTransport<DGadvection>>(
        *this->smesh, *this->meshData, *cG2DGAdvectInterpolator);
}

/*************************************************************/
template <int DGadvection>
ModelArray KokkosCGDynamicsKernel<DGadvection>::getDG0Data(const std::string& name) const
{
    HField data(ModelArray::Type::H);
    if (name == shearName) {
        computeShearDevice(tempDataAdvectDevice, e11Device, e12Device, e22Device, PSIStressDevice,
            meshData->landMaskDevice, iMJwPSIAdvectDevice);
        Kokkos::deep_copy(this->tempDataAdvectHost, this->tempDataAdvectDevice);
        return DGModelArray::dg2ma(tempDataAdvect, data);
    } else if (name == divergenceName) {
        computeTensorInvariantIDevice(tempDataDevice, e11Device, e12Device, e22Device);
        Kokkos::deep_copy(this->tempDataHost, this->tempDataDevice);
        return DGModelArray::dg2ma(tempData, data);
    } else if (name == sigmaIName) {
        computeTensorInvariantIDevice(tempDataDevice, s11Device, s12Device, s22Device);
        Kokkos::deep_copy(this->tempDataHost, this->tempDataDevice);
        return DGModelArray::dg2ma(tempData, data);
    } else if (name == sigmaIIName) {
        computeTensorInvariantIIDevice(tempDataDevice, s11Device, s12Device, s22Device);
        Kokkos::deep_copy(this->tempDataHost, this->tempDataDevice);
        return DGModelArray::dg2ma(tempData, data);
    } else {
        return CGDynamicsKernel<DGadvection>::getDG0Data(name);
    }
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::advectAndLimit(const FloatType dt,
    const ConstKokkosDeviceView<CGVector<CGdegree>>& cgUDevice,
    const ConstKokkosDeviceView<CGVector<CGdegree>>& cgVDevice)
{
    dGTransportDevice->prepareAdvection(cgUDevice, cgVDevice);

    //! Perform transport step
    dGTransportDevice->step(dt, ciceDevice);
    dGTransportDevice->step(dt, hiceDevice);

    //! Gauss-point limiting
    limitMax(ciceDevice, 1.0);
    limitMin(ciceDevice, 0.0);
    limitMin(hiceDevice, 0.0);
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::updateGradientOfSeaSurfaceHeight()
{
    // Reinit the gradient of the sea surface height. Not done by DataMap as seaSurfaceHeight is
    // always dG(0). Currently done on CPU because their are no dependencies on other
    // computations and the cost is small (<3% of dynamics with a single thread).
    this->computeGradientOfSeaSurfaceHeight(this->seaSurfaceHeight);
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(
        execSpace, this->xGradSeaSurfaceHeightDevice, this->xGradSeaSurfaceHeightHost);
    Kokkos::deep_copy(
        execSpace, this->yGradSeaSurfaceHeightDevice, this->yGradSeaSurfaceHeightHost);
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::prepareIterationDevice(const DeviceViewCG& cgHDevice,
    const DeviceViewCG& cgADevice, const ConstDeviceViewAdvect& hiceDevice,
    const ConstDeviceViewAdvect& ciceDevice,
    const Interpolations::KokkosDG2CGInterpolator<CGdegree, DGadvection>& dG2CGInterpolator)
{
    // interpolate ice height and concentration to local cg Variables
    dG2CGInterpolator(cgHDevice, hiceDevice);
    // VectorManipulations::CGAveragePeriodic(*smesh, cgH);
    dG2CGInterpolator(cgADevice, ciceDevice);
    // VectorManipulations::CGAveragePeriodic(*smesh, cgA);

    /* limit A to [0,1] and H to [5 cm, ...)
     * This limit on H is equivalent to assuming that ice thinner than 5 cm is always in free
     * drift, which is reasonable. We need a limit of the order of cm here, so that the solver
     * remains stable. With a limit of the order of mm, we need a much smaller time step to
     * remain stable.
     */
    Kokkos::parallel_for(
        "limitCGA", cgADevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex idx) {
            cgADevice(idx) = Kokkos::fmin(Kokkos::fmax(cgADevice(idx), static_cast<FloatType>(0.0)),
                static_cast<FloatType>(1.0));
        });
    Kokkos::parallel_for(
        "limitCGH", cgHDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex idx) {
            cgHDevice(idx) = Kokkos::fmax(cgHDevice(idx), static_cast<FloatType>(0.05));
        });
}

/*************************************************************/
template <int CG, typename Vec>
static KOKKOS_IMPL_FUNCTION Eigen::Matrix<FloatType, CGDOFS(CG), 1> cgToLocal(
    const Vec& vGlobal, DeviceIndex cgi, DeviceIndex cgShift)
{
    if constexpr (CG == 1) {
        Eigen::Matrix<FloatType, CGDOFS(1), 1> vLocal;
        vLocal << vGlobal(cgi), vGlobal(cgi + 1), vGlobal(cgi + cgShift),
            vGlobal(cgi + 1 + cgShift);
        return vLocal;
    } else {
        Eigen::Matrix<FloatType, CGDOFS(2), 1> vLocal;
        vLocal << vGlobal(cgi), vGlobal(cgi + 1), vGlobal(cgi + 2), vGlobal(cgi + cgShift),
            vGlobal(cgi + 1 + cgShift), vGlobal(cgi + 2 + cgShift), vGlobal(cgi + 2 * cgShift),
            vGlobal(cgi + 1 + 2 * cgShift), vGlobal(cgi + 2 + 2 * cgShift);
        return vLocal;
    }
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::projectVelocityToStrainDevice(
    const ConstDeviceViewCG& uDevice, const ConstDeviceViewCG& vDevice,
    const DeviceViewStress& e11Device, const DeviceViewStress& e12Device,
    const DeviceViewStress& e22Device, const ConstDeviceBitset& landMaskDevice,
    const GradMapDevice& iMgradXDevice, const GradMapDevice& iMgradYDevice,
    const GradMapDevice& iMMDevice, DeviceIndex nx, DeviceIndex ny, COORDINATES coordinates)
{
    const DeviceIndex cgshift = CGdegree * nx + 1; //!< Index shift for each row

    // 1D loop is much faster than 2D loop on CPU
    Kokkos::parallel_for(
        "projectVelocityToStrain", nx * ny, KOKKOS_LAMBDA(const DeviceIndex idx) {
            // only on ice
            if (!landMaskDevice.test(idx)) {
                return;
            }

            const DeviceIndex col = idx % nx;
            const DeviceIndex row = idx / nx;
            const DeviceIndex dgi = nx * row + col; //!< Index of dg vector
            const DeviceIndex cgi
                = CGdegree * cgshift * row + col * CGdegree; //!< Lower left index of cg vector

            const auto u = makeEigenMap(uDevice);
            const auto v = makeEigenMap(vDevice);

            auto e11 = makeEigenMap(e11Device);
            auto e12 = makeEigenMap(e12Device);
            auto e22 = makeEigenMap(e22Device);

            // get the local x/y - velocity coefficients on the element
            const auto vxLocal = cgToLocal<CGdegree>(u, cgi, cgshift);
            const auto vyLocal = cgToLocal<CGdegree>(v, cgi, cgshift);

            // Solve (E, Psi) = (0.5(DV + DV^T), Psi)
            // by integrating rhs and inverting with dG(stress) mass matrix
            const auto iMgradX = iMgradXDevice[dgi];
            e11.row(dgi) = iMgradX * vxLocal;
            const auto iMgradY = iMgradYDevice[dgi];
            e22.row(dgi) = iMgradY * vyLocal;
            e12.row(dgi) = 0.5 * (iMgradX * vyLocal + iMgradY * vxLocal);

            if (coordinates == SPHERICAL) {
                const auto iMM = iMMDevice[dgi];
                e11.row(dgi) -= iMM * vyLocal;
                e12.row(dgi) += 0.5 * iMM * vxLocal;
            }
        });
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::dirichletZero(
    const DeviceViewCG& v, const ConstDeviceBitset& cgLandMaskDevice)
{
    // TR 07.04.2025: Dirichlet Zero (u=v=0) holds on land and on the boundary between
    // land and ice. Hence on all elements with landmask = 0, or, on cg nodes with
    // cglandmask = 0
    Kokkos::parallel_for(
        "dirichletZero", v.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            if (!cgLandMaskDevice.test(i)) {
                v(i) = 0.0;
            }
        });
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::computeStressDivergenceDevice(
    const DeviceViewCG& dStressXDevice, const DeviceViewCG& dStressYDevice,
    const ConstDeviceViewStress& s11Device, const ConstDeviceViewStress& s12Device,
    const ConstDeviceViewStress& s22Device, const ConstDeviceBitset& landMaskDevice,
    const DivMapDevice& divS1Device, const DivMapDevice& divS2Device,
    const DivMapDevice& divMDevice, const ConstDeviceBitset& cgLandMaskDevice, DeviceIndex nx,
    DeviceIndex ny, COORDINATES coordinates)
{
    using CGVec = Eigen::Vector<Nextsim::FloatType, CGdof>;
    //    static PerfTimer timerDivZero("divZeroGPU");
    // static PerfTimer timerDivComp("divCompGPU");
    //   static PerfTimer timerDivDirichlet("divDirichletGPU");
    // zero buffers
    //   timerDivZero.start();
    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, dStressXDevice, 0.0);
    Kokkos::deep_copy(execSpace, dStressYDevice, 0.0);
    //  execSpace.fence();
    //    timerDivZero.stop();

    //   timerDivComp.start();
    Kokkos::parallel_for(
        "computeStressDivergence", nx * ny, KOKKOS_LAMBDA(const DeviceIndex idx) {
            // only on ice!
            if (!landMaskDevice.test(idx)) {
                return;
            }

            const DeviceIndex cx = idx % nx;
            const DeviceIndex cy = idx / nx;
            const DeviceIndex eid = cx + nx * cy;

            const auto s11 = makeEigenMap(s11Device);
            const auto s12 = makeEigenMap(s12Device);
            const auto s22 = makeEigenMap(s22Device);

            const auto divS1 = divS1Device[eid];
            const auto divS2 = divS2Device[eid];
            CGVec tx = divS1 * s11.row(eid).transpose() + divS2 * s12.row(eid).transpose();
            CGVec ty = divS1 * s12.row(eid).transpose() + divS2 * s22.row(eid).transpose();

            if (coordinates == SPHERICAL) {
                const auto divM = divMDevice[eid];
                tx += divM * s12.row(eid).transpose();
                ty -= divM * s11.row(eid).transpose();
            }

            const DeviceIndex cgRow = CGdegree * nx + 1;
            //!< lower left CG-index in element (cx,cy)
            const DeviceIndex cg_i = CGdegree * cgRow * cy + CGdegree * cx;

            // Fill the stress divergence values
            for (DeviceIndex row = 0; row <= CGdegree; ++row) {
                for (DeviceIndex col = 0; col <= CGdegree; ++col) {
                    const DeviceIndex dst_idx = cg_i + col + row * cgRow;
                    const DeviceIndex src_idx = col + (CGdegree + 1) * row;
                    Kokkos::atomic_sub(&dStressXDevice(dst_idx), tx(src_idx));
                    Kokkos::atomic_sub(&dStressYDevice(dst_idx), ty(src_idx));
                }
            }
        });
    //    timerDivComp.stop();
    // set zero on the Dirichlet boundaries
    //   timerDivDirichlet.start();
    dirichletZero(dStressXDevice, cgLandMaskDevice);
    dirichletZero(dStressYDevice, cgLandMaskDevice);
    //  timerDivDirichlet.stop();

    // todo: add the contributions on the periodic boundaries
    // VectorManipulations::CGAveragePeriodic(*smesh, tx);
    //   VectorManipulations::CGAveragePeriodic(*smesh, ty);
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::applyBoundariesDevice(const DeviceViewCG& uDevice,
    const DeviceViewCG& vDevice, const ConstDeviceBitset& cgLandMaskDevice, DeviceIndex nx,
    DeviceIndex ny)
{
    dirichletZero(uDevice, cgLandMaskDevice);
    dirichletZero(vDevice, cgLandMaskDevice);

    // TODO Periodic boundary conditions.
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::updateIceOceanStressDevice(
    const DeviceViewCG& uIceOceanStressDevice, const DeviceViewCG& vIceOceanStressDevice,
    const ConstDeviceViewCG& uIceDevice, const ConstDeviceViewCG& vIceDevice,
    const ConstDeviceViewCG& uOceanDevice, const ConstDeviceViewCG& vOceanDevice,
    const DynamicsParameters& params, FloatType cosOceanAngle, FloatType sinOceanAngle)
{
    const FloatType FOcean = params.COcean * params.rhoOcean;

    Kokkos::parallel_for(
        "computeStressDivergence", uIceOceanStressDevice.extent(0),
        KOKKOS_LAMBDA(const DeviceIndex i) {
            const FloatType uOceanRel = uOceanDevice(i) - uIceDevice(i);
            const FloatType vOceanRel = vOceanDevice(i) - vIceDevice(i);
            const FloatType cPrime = FOcean * Kokkos::hypot(uOceanRel, vOceanRel);
            uIceOceanStressDevice(i)
                = cPrime * (uOceanRel * cosOceanAngle - vOceanRel * sinOceanAngle);
            vIceOceanStressDevice(i)
                = cPrime * (vOceanRel * cosOceanAngle + uOceanRel * sinOceanAngle);
        });
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::computeShearDevice(const DeviceViewAdvect& destDevice,
    const ConstDeviceViewStress& e11Device, const ConstDeviceViewStress& e12Device,
    const ConstDeviceViewStress& e22Device, const PSIStressView& PSIStressDevice,
    const ConstDeviceBitset& landMaskDevice, const GaussMapAdvectDevice& iMJwPSIAdvectDevice)
{
    Kokkos::parallel_for(
        "computeShearDevice", destDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto dest = makeEigenMap(destDevice);
            if (!landMaskDevice.test(i)) {
                dest.row(i).setZero();
                return;
            }

            const auto e11 = makeEigenMap(e11Device);
            const auto e12 = makeEigenMap(e12Device);
            const auto e22 = makeEigenMap(e22Device);

            const auto PSIStress = makeEigenMap(PSIStressDevice);

            const EdgeVec e11Gauss = e11.row(i) * PSIStress;
            const EdgeVec e12Gauss = e12.row(i) * PSIStress;
            const EdgeVec e22Gauss = e22.row(i) * PSIStress;

            dest.row(i) = iMJwPSIAdvectDevice[i]
                * ((e11Gauss.array() - e22Gauss.array()).square() + 4.0 * e12Gauss.array().square()
                    + 1.e-20)
                      .sqrt()
                      .matrix()
                      .transpose();
        });
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::computeTensorInvariantIDevice(
    const DeviceViewStress& destDevice, const ConstDeviceViewStress& e11Device,
    const ConstDeviceViewStress& e12Device, const ConstDeviceViewStress& e22Device)
{
    Kokkos::parallel_for(
        "computeTensorInvariant", destDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto dest = makeEigenMap(destDevice);
            const auto e11 = makeEigenMap(e11Device);
            const auto e12 = makeEigenMap(e12Device);
            const auto e22 = makeEigenMap(e22Device);

            dest.row(i) = 0.5 * (e11.row(i) + e22.row(i));
        });
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::computeTensorInvariantIIDevice(
    const DeviceViewStress& destDevice, const ConstDeviceViewStress& e11Device,
    const ConstDeviceViewStress& e12Device, const ConstDeviceViewStress& e22Device)
{
    Kokkos::parallel_for(
        "computeTensorInvariant", destDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto dest = makeEigenMap(destDevice);
            const auto e11 = makeEigenMap(e11Device);
            const auto e12 = makeEigenMap(e12Device);
            const auto e22 = makeEigenMap(e22Device);

            dest.row(i)
                = ((0.5 * (e11.row(i) - e22.row(i))).array().square() + e12.row(i).array().square())
                      .sqrt();
        });
}

/*************************************************************/
// because ParametricMomentumMap<CGdegree>::iMJwPSIAdvect does not properly depend on DGadvection we
// can only build this version since the switch is implemented in compile-time we dont really need
// the other versions anyway
template class KokkosCGDynamicsKernel<DGCOMP>;
// template class KokkosCGDynamicsKernel<1>;
// template class KokkosCGDynamicsKernel<3>;
// template class KokkosCGDynamicsKernel<6>;
}