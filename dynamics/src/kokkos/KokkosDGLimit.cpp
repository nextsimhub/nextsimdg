/*!
 * @file dgLimit.cpp
 * @date 27 September 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosDGLimit.hpp"
#include "../include/codeGenerationDGinGauss.hpp"

namespace Nextsim {

void limitMax(const KokkosDeviceView<DGVector<1>>& dg, FloatType max)
{
    Kokkos::parallel_for(
        "limitMax", dg.size(),
        KOKKOS_LAMBDA(const DeviceIndex i) { dg(i) = std::min(dg(i), max); });
}

void limitMax(const KokkosDeviceView<DGVector<3>>& dg, FloatType max)
{
    Kokkos::parallel_for(
        "limitMax", dg.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            dg(i, 0) = std::min(max, dg(i, 0));
            const FloatType l0 = static_cast<FloatType>(2.0)
                * std::max(Kokkos::abs(dg(i, 1) + dg(i, 2)), Kokkos::abs(dg(i, 1) - dg(i, 2)));
            if (l0 == 0)
                return;
            const FloatType ex = dg(i, 0) + l0 - max;
            if (ex > 0) {
                dg(i, 1) *= (max - dg(i, 0)) / l0;
                dg(i, 2) *= (max - dg(i, 0)) / l0;
            }
        });
}

void limitMax(const KokkosDeviceView<DGVector<6>>& dgDevice, FloatType max)
{
    const auto PSI63 = PSI<6, 3>;
    Kokkos::parallel_for(
        "limitMax", dgDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            auto dg = makeEigenMap(dgDevice);

            dg(i, 0) = std::min(max, dg(i, 0));

            const FloatType maxvalue = (dg.block<1, 6>(i, 0) * PSI63).maxCoeff();

            // part coming from the nonlinearity
            const FloatType l0 = maxvalue - dg(i, 0);
            // value to big?
            if (maxvalue > max) {
                // limit nonlinear part
                dg.block<1, 5>(i, 1) *= (max - dg(i, 0)) / l0;
            }
        });
}

/*************************************************************/
void limitMin(const KokkosDeviceView<DGVector<1>>& dg, FloatType min)
{
    Kokkos::parallel_for(
        "limitMin", dg.size(),
        KOKKOS_LAMBDA(const DeviceIndex i) { dg(i) = std::max(dg(i), min); });
}

void limitMin(const KokkosDeviceView<DGVector<3>>& dg, FloatType min)
{
    Kokkos::parallel_for(
        "limitMin", dg.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
            dg(i, 0) = std::max(min, dg(i, 0));
            const FloatType l0 = static_cast<FloatType>(2.0)
                * std::max(Kokkos::abs(dg(i, 1) + dg(i, 2)), Kokkos::abs(dg(i, 1) - dg(i, 2)));
            if (l0 == 0)
                return;

            const FloatType ex = dg(i, 0) - l0 - min;
            if (ex < 0) {
                dg(i, 1) *= (dg(i, 0) - min) / l0;
                dg(i, 2) *= (dg(i, 0) - min) / l0;
            }
        });
}
void limitMin(const KokkosDeviceView<DGVector<6>>& dgDevice, FloatType min)
{
	const auto PSI63 = PSI<6, 3>;
    Kokkos::parallel_for(
        "limitMin", dgDevice.extent(0), KOKKOS_LAMBDA(const DeviceIndex i) {
			auto dg = makeEigenMap(dgDevice);

            dg(i, 0) = std::max(min, dg(i, 0));

            const FloatType minvalue = (dg.block<1, 6>(i, 0) * PSI63).minCoeff();

            // part coming from nonlinearity
            const FloatType l0 = minvalue - dg(i, 0);
            if (minvalue < min) {
                dg.block<1, 5>(i, 1) *= (dg(i, 0) - min) / l0;
            }
        });
}
}