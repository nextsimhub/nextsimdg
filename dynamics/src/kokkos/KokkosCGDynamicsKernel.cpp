/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosCGDynamicsKernel.hpp"
#include "../../../core/src/kokkos/include/KokkosTimer.hpp"
#include "include/KokkosDGLimit.hpp"
#include "include/KokkosDGTransport.hpp"
#include "include/KokkosInterpolations.hpp"
#include "include/KokkosMesh.hpp"
#include "include/KokkosSlopeLimiter.hpp"

namespace Nextsim {
/*************************************************************/
template <int DGadvection>
KokkosCGDynamicsKernel<DGadvection>::KokkosCGDynamicsKernel(const DynamicsParameters& params)
    : CGDynamicsKernel<DGadvection>(params)
{
}

// defined explicitly to enable pimpl with unique ptr
template <int DGadvection> KokkosCGDynamicsKernel<DGadvection>::~KokkosCGDynamicsKernel() = default;

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
    std::tie(seaSurfaceHeightHost, seaSurfaceHeightDevice)
        = makeKokkosDualView("seaSurfaceHeight", this->seaSurfaceHeight);

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
    ModelArray tempDataMA;
    tempDataMADevice = makeKokkosDeviceView("tempDataMA", tempDataMA.data());

    assert(this->pmap);
    divS1Device = makeKokkosDeviceViewMap("divS1", this->pmap->divS1, MakeViewOptions::DEVICE_COPY);
    divS2Device = makeKokkosDeviceViewMap("divS2", this->pmap->divS2, MakeViewOptions::DEVICE_COPY);
    divMDevice = makeKokkosDeviceViewMap("divM", this->pmap->divM, MakeViewOptions::DEVICE_COPY);
    iMgradXDevice
        = makeKokkosDeviceViewMap("iMgradX", this->pmap->iMgradX, MakeViewOptions::DEVICE_COPY);
    iMgradYDevice
        = makeKokkosDeviceViewMap("iMgradY", this->pmap->iMgradY, MakeViewOptions::DEVICE_COPY);
    iMMDevice = makeKokkosDeviceViewMap("iMM", this->pmap->iMM, MakeViewOptions::DEVICE_COPY);

    PSIAdvectDevice = makeKokkosDeviceView(
        "PSI<DGadvection, NGP>", PSI<DGadvection, NGP>, MakeViewOptions::DEVICE_COPY);
    PSIStressDevice = makeKokkosDeviceView(
        "PSI<DGstress, NGP>", PSI<DGstressComp, NGP>, MakeViewOptions::DEVICE_COPY);

    lumpedCGMassDevice = makeKokkosDeviceView(
        "lumpedCGMass", this->pmap->lumpedcgmass, MakeViewOptions::DEVICE_COPY);
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
    slopeLimiterDevice
        = std::make_unique<KokkosSlopeLimiter<DGadvection>>(*this->smesh, *this->meshData);

    if constexpr (DGadvection != 1 || CGdegree != 1) {
        dG2CGFirstOrderInterpolator
            = std::make_unique<Interpolations::KokkosDG2CGInterpolator<1, 1>>(*this->smesh);
        _dXSSHDevice
            = makeKokkosDeviceViewMap("dXSSH", this->pmap->dX_SSH, MakeViewOptions::DEVICE_COPY);
        _dYSSHDevice
            = makeKokkosDeviceViewMap("dYSSH", this->pmap->dY_SSH, MakeViewOptions::DEVICE_COPY);
        _lumpedCG1MassDevice = makeKokkosDeviceView(
            "lumpedCG1Mass", this->pmap->lumpedcg1mass, MakeViewOptions::DEVICE_COPY);
        _uGradDevice = DeviceViewCG1("uGrad", _lumpedCG1MassDevice.extent(0));
        _vGradDevice = DeviceViewCG1("vGrad", _lumpedCG1MassDevice.extent(0));
    }

    namedCGFields = {
        { uName, uDevice },
        { vName, vDevice },
        { uWindName, uAtmosDevice },
        { vWindName, vAtmosDevice },
        { uOceanName, uOceanDevice },
        { vOceanName, vOceanDevice },
    };
}

/*************************************************************/
template <int DGadvection>
ModelArray KokkosCGDynamicsKernel<DGadvection>::getDG0Data(const std::string& name) const
{
    if (name == shearName) {
        computeShearDevice(tempDataAdvectDevice, e11Device, e12Device, e22Device, PSIStressDevice,
            meshData->landMaskDevice, iMJwPSIAdvectDevice);
        HField data(ModelArray::Type::H);
        return dG2MA(data, tempDataAdvectDevice);
    } else if (name == divergenceName) {
        computeTensorInvariantIDevice(tempDataDevice, e11Device, e12Device, e22Device);
        Kokkos::deep_copy(this->tempDataHost, this->tempDataDevice);
        HField data(ModelArray::Type::H);
        return dG2MA(data, tempDataAdvectDevice);
    } else if (name == sigmaIName) {
        computeTensorInvariantIDevice(tempDataDevice, s11Device, s12Device, s22Device);
        HField data(ModelArray::Type::H);
        return dG2MA(data, tempDataAdvectDevice);
    } else if (name == sigmaIIName) {
        computeTensorInvariantIIDevice(tempDataDevice, s11Device, s12Device, s22Device);
        HField data(ModelArray::Type::H);
        return dG2MA(data, tempDataAdvectDevice);
    } else if (name == uName) {
        (*cG2DGAdvectInterpolator)(tempDataAdvectDevice, uDevice);
        HField data(ModelArray::Type::U);
        return dG2MA(data, tempDataAdvectDevice);
    } else if (name == vName) {
        (*cG2DGAdvectInterpolator)(tempDataAdvectDevice, vDevice);
        HField data(ModelArray::Type::V);
        return dG2MA(data, tempDataAdvectDevice);
    } else if (name == uIOStressName) {
        (*cG2DGAdvectInterpolator)(tempDataAdvectDevice, uIceOceanStressDevice);
        HField data(ModelArray::Type::U);
        return dG2MA(data, tempDataAdvectDevice);
    } else if (name == vIOStressName) {
        (*cG2DGAdvectInterpolator)(tempDataAdvectDevice, vIceOceanStressDevice);
        HField data(ModelArray::Type::V);
        return dG2MA(data, tempDataAdvectDevice);
    } else {
        return CGDynamicsKernel<DGadvection>::getDG0Data(name);
    }
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::setData(const std::string& name, const ModelArray& data)
{
    // Just copy to device and use the other setter so that we don't have to treat every field
    // differently. This is much faster than directly copying to the destination if there is
    // a stride involved.
    // const auto& [dataHost, dataDevice] = makeKokkosDualView(name + "Temp", data.data());
    assert(data.components(0).size() == 1 && "Expecting only dg(0) fields in setData()");
    const auto dataHost = makeKokkosHostView(data.data());
    Kokkos::deep_copy(tempDataMADevice, dataHost);
    setData(name, tempDataMADevice);
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::setData(
    const std::string& name, const ConstDeviceViewMA& data)
{
    // Special cases: hice, cice
    if (name == hiceName || name == ciceName || name == hsnowName) {
        throw std::runtime_error(std::string("Use setDGArray() to set the data for ") + name);
    } else if (name == sshName) {
        kokkosMA2DG<1>(seaSurfaceHeightDevice, data);
    } else if (auto it = namedCGFields.find(name); it != namedCGFields.end()) {
        mA2CG(it->second, data);
    } else {
        throw std::runtime_error(std::string("Trying to setData() for the unknown field ") + name);
    }
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::setDGArray(
    const std::string& name, ModelArray::DataType& dgData)
{
    // There are different problems depending on the configuration.
    if constexpr (IS_GPU_EXEC_SPACE<Kokkos::DefaultExecutionSpace>) {
        // If the kernel runs on the device it would be necessary to copy every DGVector set with
        // this function back after the update. It makes more sense to handle host-device transfers
        // outside of the kernel where the ModelArrayStore tracks the state.
        throw std::runtime_error("Setting the buffers with setDGArray(ModelArray::DataType&) does "
                                 "not work properly with device execution "
                                 "because the updated fields are not transferred back.");
    } else {
        // To implement this function we just have to update the correct Kokkos views, i.e.
        // std::tie(hiceHost, hiceDevice) = makeKokkosDualView("hice", dgData);
        // However, there is currently no scenario where this function is needed.
        // As long as Kokkos is enabled, the view based setDGArray() will be used, even without
        // device execution.
        throw std::runtime_error(
            "Setting the buffers with setDGArray(ModelArray::DataType&) is currently"
            "not implemented for the Kokkos kernel.");
    }
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::setDGArray(
    const std::string& name, const KokkosDeviceView<ModelArray::DataType>& dgData)
{
    if (name == hiceName) {
        hiceDevice = dgData;
    } else if (name == ciceName) {
        ciceDevice = dgData;
    } else if (name == hsnowName) {
        hsnowDevice = dgData;
    }
}

/*************************************************************/
template <int DGadvection> void KokkosCGDynamicsKernel<DGadvection>::prepareAdvection()
{
    static KokkosTimer<DETAILED_MEASUREMENTS> timerPrepAdvection("prepareAdvection");
    timerPrepAdvection.start();
    dGTransportDevice->prepareAdvection(uDevice, vDevice);
    timerPrepAdvection.stop();
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::advectDynamicsFields(FloatType timestep)
{
    advectDGVFieldDevice(timestep, hiceDevice, 0.0);
    advectDGVFieldDevice(timestep, ciceDevice, 0.0, 1.0);
    advectDGVFieldDevice(timestep, hsnowDevice, 0.0);
}

/*************************************************************/
template <int DGadvection>
DGVector<DGadvection>& KokkosCGDynamicsKernel<DGadvection>::advectDGVField(
    FloatType timestep, DGVector<DGadvection>& field, FloatType lowerLimit, FloatType upperLimit)
{
    static KokkosTimer<DETAILED_MEASUREMENTS> timer("advectExternalGPU");

    timer.start();
    auto fieldViewHost = makeKokkosHostView(field);
    Kokkos::deep_copy(this->tempDataAdvectDevice, fieldViewHost);
    advectDGVFieldDevice(timestep, this->tempDataAdvectDevice, lowerLimit, upperLimit);
    Kokkos::deep_copy(fieldViewHost, this->tempDataAdvectDevice);
    timer.stop();

    return field;
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::advectDGVFieldDevice(
    FloatType timestep, const DeviceViewAdvect& field, FloatType lowerLimit, FloatType upperLimit)
{
    dGTransportDevice->step(timestep, field);

    //! Slope Limiting
    bool limitSlope = false;

    // First, limit minimum and/or maximum of the average component
    if (lowerLimit > -std::numeric_limits<FloatType>::infinity()) {
        slopeLimiterDevice->limitMin(field, lowerLimit);
        limitSlope = true;
    }
    if (upperLimit < std::numeric_limits<FloatType>::infinity()) {
        slopeLimiterDevice->limitMax(field, upperLimit);
        limitSlope = true;
    }

    // Then prevent new local minima and maxima
    if (limitSlope)
        slopeLimiterDevice->limit(field);
}

/*************************************************************/
template <typename Mat> void compare(const std::string& name, const Mat& m1, const Mat& m2)
{
    FloatType normRef = m1.norm();
    FloatType normDiff = (m1 - m2).norm();
    std::cout << name << " - abs: " << normDiff << ", rel: " << normDiff / normRef
              << ", norm: " << normRef;
    Eigen::Index maxIndex;
    const FloatType maxVal = (m1 - m2).cwiseAbs().maxCoeff(&maxIndex);
    std::cout << ", max diff: " << maxVal << " at " << maxIndex << " abs(" << m1(maxIndex) << "-"
              << m2(maxIndex) << ")" << std::endl;
}

template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::updateGradientOfSeaSurfaceHeight()
{
    // capturing structured bindings is C++20 so we can't just do
    // const auto& [dG2CGInterpolater, dXSSHDevice, dYSSHDevice, lumpedCG1MassDevice] =
    const auto precomputedMaps = [this]() {
        if constexpr (DGadvection == 1 && CGdegree == 1) {
            return std::tie(dG2CGAdvectInterpolator, divS1Device, divS2Device, lumpedCGMassDevice);
        } else {
            return std::tie(
                dG2CGFirstOrderInterpolator, _dXSSHDevice, _dYSSHDevice, _lumpedCG1MassDevice);
        }
    }();
    const auto& dG2CGInterpolater = std::get<0>(precomputedMaps);
    const auto& dXSSHDevice = std::get<1>(precomputedMaps);
    const auto& dYSSHDevice = std::get<2>(precomputedMaps);
    const auto& lumpedCG1MassDevice = std::get<3>(precomputedMaps);

    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, this->seaSurfaceHeightDevice, this->seaSurfaceHeightHost);

    const DeviceViewCG1 cgSeaSurfaceHeightDevice
        = DeviceViewCG1("cgSeaSurfaceHeight", lumpedCG1MassDevice.extent(0));
    (*dG2CGInterpolater)(cgSeaSurfaceHeightDevice, seaSurfaceHeightDevice);

    const auto gradientFields = [&]() {
        if constexpr (CGdegree == 1) {
            return std::tie(xGradSeaSurfaceHeightDevice, yGradSeaSurfaceHeightDevice);
        } else {
            return std::tie(_uGradDevice, _vGradDevice);
        }
    }();
    const DeviceViewCG1& uGradDevice = std::get<0>(gradientFields);
    const DeviceViewCG1& vGradDevice = std::get<1>(gradientFields);
    Kokkos::deep_copy(execSpace, uGradDevice, 0.0);
    Kokkos::deep_copy(execSpace, vGradDevice, 0.0);

    const DeviceIndex nx = this->smesh->nx;
    const DeviceIndex ny = this->smesh->ny;

    Kokkos::parallel_for(
        "computeCGSeaSurfaceGrads", nx * ny, KOKKOS_LAMBDA(const DeviceIndex eid) {
            const DeviceIndex cy = eid / nx; //!< y-index of element
            const DeviceIndex cx = eid % nx; //!< x-index of element
            const DeviceIndex cg1id = cy * (nx + 1) + cx; //!< lower/left Index in cg vector

            // get local CG nodes
            const Eigen::Vector<FloatType, 4> locCGSSH = { cgSeaSurfaceHeightDevice(cg1id),
                cgSeaSurfaceHeightDevice(cg1id + 1), cgSeaSurfaceHeightDevice(cg1id + nx + 1),
                cgSeaSurfaceHeightDevice(cg1id + nx + 1 + 1) };

            const auto dXSSH = makeEigenMap(dXSSHDevice);
            const auto dYSSH = makeEigenMap(dYSSHDevice);

            // compute grad
            const Eigen::Vector<FloatType, 4> tx = dXSSH[eid] * locCGSSH;
            const Eigen::Vector<FloatType, 4> ty = dYSSH[eid] * locCGSSH;

            // add global vector
            Kokkos::atomic_sub(&uGradDevice(cg1id), tx(0));
            Kokkos::atomic_sub(&uGradDevice(cg1id + 1), tx(1));
            Kokkos::atomic_sub(&uGradDevice(cg1id + nx + 1), tx(2));
            Kokkos::atomic_sub(&uGradDevice(cg1id + nx + 1 + 1), tx(3));
            Kokkos::atomic_sub(&vGradDevice(cg1id), ty(0));
            Kokkos::atomic_sub(&vGradDevice(cg1id + 1), ty(1));
            Kokkos::atomic_sub(&vGradDevice(cg1id + nx + 1), ty(2));
            Kokkos::atomic_sub(&vGradDevice(cg1id + nx + 1 + 1), ty(3));
        });

    // scale with mass
    const DeviceIndex cg1Row = nx + 1;
    // const DeviceIndex cg1rowRed = nx - 1;
    const auto makeScaleGradFn = [&](const DeviceViewCG1& gradDevice) {
        return KOKKOS_LAMBDA(const DeviceIndex i) { gradDevice(i) /= lumpedCG1MassDevice(i); };
    };

    Kokkos::parallel_for("scaleUGrad", uGradDevice.extent(0), makeScaleGradFn(uGradDevice));
    Kokkos::parallel_for("scaleVGrad", vGradDevice.extent(0), makeScaleGradFn(vGradDevice));

    // correct boundary (just extend in last elements)
    // Corners are handled impliclty. By treating them like regular nodes they recieve the
    // values of the inner diagonal neighbors during the second (y) update.
    const DeviceIndex topLeft = ny * cg1Row;
    const auto makeExtendBoundaryXFn = [&](const DeviceViewCG1& gradDevice) {
        return KOKKOS_LAMBDA(const DeviceIndex i)
        {
            gradDevice(i) = gradDevice(i + cg1Row);
            gradDevice(topLeft + i) = gradDevice(topLeft + i - cg1Row);
        };
    };

    const auto makeExtendBoundaryYFn = [&](const DeviceViewCG1& gradDevice) {
        return KOKKOS_LAMBDA(const DeviceIndex i)
        {
            const DeviceIndex j = i * cg1Row;
            gradDevice(j) = gradDevice(j + 1);
            gradDevice(j + cg1Row - 1) = gradDevice(j + cg1Row - 1 - 1);
        };
    };

    Kokkos::parallel_for("extendCornersUGradX", nx + 1, makeExtendBoundaryXFn(uGradDevice));
    Kokkos::parallel_for("extendCornersUGradY", ny + 1, makeExtendBoundaryYFn(uGradDevice));

    Kokkos::parallel_for("extendCornersVGradX", nx + 1, makeExtendBoundaryXFn(vGradDevice));
    Kokkos::parallel_for("extendCornersVGradY", ny + 1, makeExtendBoundaryYFn(vGradDevice));

    // Interpolate to CG2
    // If we have CGdegree == 1 we are already finished since the computations where done
    // directly in the destination (xGradSeaSurfaceHeightDevice, yGradSeaSurfaceHeightDevice).
    if constexpr (CGdegree == 2) {
        Kokkos::deep_copy(execSpace, xGradSeaSurfaceHeightDevice, 0.0);
        Kokkos::deep_copy(execSpace, yGradSeaSurfaceHeightDevice, 0.0);
        Interpolations::kokkosCG12CG2(xGradSeaSurfaceHeightDevice, uGradDevice, nx, ny);
        Interpolations::kokkosCG12CG2(yGradSeaSurfaceHeightDevice, vGradDevice, nx, ny);
    }

    // gradients are only used internally so they don't have to be synced
    /*    Kokkos::deep_copy(
            execSpace, this->xGradSeaSurfaceHeightHost, xGradSeaSurfaceHeightDevice);
        Kokkos::deep_copy(
            execSpace, this->yGradSeaSurfaceHeightHost, yGradSeaSurfaceHeightDevice);*/
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
            e12.row(dgi) = FloatType(0.5) * (iMgradX * vyLocal + iMgradY * vxLocal);

            if (coordinates == SPHERICAL) {
                const auto iMM = iMMDevice[dgi];
                e11.row(dgi) -= iMM * vyLocal;
                e12.row(dgi) += FloatType(0.5) * iMM * vxLocal;
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
    using CGVec = Eigen::Vector<FloatType, CGdof>;

    auto execSpace = Kokkos::DefaultExecutionSpace();
    Kokkos::deep_copy(execSpace, dStressXDevice, 0.0);
    Kokkos::deep_copy(execSpace, dStressYDevice, 0.0);

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
    // set zero on the Dirichlet boundaries
    dirichletZero(dStressXDevice, cgLandMaskDevice);
    dirichletZero(dStressYDevice, cgLandMaskDevice);
}

/*************************************************************/
template <int DGadvection>
void KokkosCGDynamicsKernel<DGadvection>::applyBoundariesDevice(const DeviceViewCG& uDevice,
    const DeviceViewCG& vDevice, const ConstDeviceBitset& cgLandMaskDevice, DeviceIndex nx,
    DeviceIndex ny)
{
    dirichletZero(uDevice, cgLandMaskDevice);
    dirichletZero(vDevice, cgLandMaskDevice);
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

            dest.row(i) = FloatType(0.5) * (e11.row(i) + e22.row(i));
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

            dest.row(i) = ((FloatType(0.5) * (e11.row(i) - e22.row(i))).array().square()
                + e12.row(i).array().square())
                              .sqrt();
        });
}

/*************************************************************/
template <int DGadvection>
ModelArray& KokkosCGDynamicsKernel<DGadvection>::dG2MA(
    ModelArray& ma, const ConstDeviceViewAdvect& dg) const
{
    kokkosDG2MA<DGadvection>(tempDataMADevice, dg);
    const auto maView = makeKokkosHostView(ma.getDataRef());
    Kokkos::deep_copy(maView, tempDataMADevice);

    return ma;
}

/*************************************************************/
// because ParametricMomentumMap<CGdegree>::iMJwPSIAdvect does not properly depend on DGadvection we
// can only build this version since the switch is implemented in compile-time we dont really need
// the other versions anyway
template class KokkosCGDynamicsKernel<DGCOMP>;

}
