/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

// by also guarding for USE_KOKKOS this header can be safely included even
// when Kokkos is not enabled
#if !defined(KOKKOSMODELARRAY_HPP) && defined(USE_KOKKOS)
#define KOKKOSMODELARRAY_HPP

#include "../../include/ModelArray.hpp"
#include "KokkosTimer.hpp"
#include "KokkosUtils.hpp"
#include "dgVector.hpp"

namespace Nextsim {

// kokkos views compatible with ModelArray
using DeviceViewMA = KokkosDeviceView<ModelArray::DataType>;
using ConstDeviceViewMA = ConstKokkosDeviceView<ModelArray::DataType>;
using HostViewMA = KokkosHostView<ModelArray::DataType>;
using ConstHostViewMA = ConstKokkosHostView<ModelArray::DataType>;

/*!
 * @brief Copy data from an ModelArray to a DGVector through Kokkos views.
 *
 * @details Source and destination do not have to be located in the same memory space but strided
 * host -> device transfers can be extremely slow.
 *
 * @param dg The destination DGVector view.
 * @param ma The source ModelArray view.
 */
template <int N, typename... ArgsMA, typename... ArgsDG>
void kokkosMA2DG(const KokkosEigenView<DGVector<N>, ArgsDG...>& dg,
    const ConstKokkosEigenView<ModelArray::DataType, ArgsMA...>& ma)
{
    // N == 1 needs different treatment at compile-time because the corresponding Kokkos::View has
    // only one dimension.
    if constexpr (N == 1) {
        assert(ma.extent(1) == 1);
        // views need to have the same rank so we squeeze 1-sized component dimension
        const auto firstCompMA = Kokkos::subview(ma, Kokkos::ALL(), 0);
        Kokkos::deep_copy(dg, firstCompMA);
    } else if (N == ma.extent(1)) {
        Kokkos::deep_copy(dg, ma);
    } else {
        assert(ma.extent(1) == 1);
        // Assign only to the 0th component. Use a range to keep the dimension.
        const auto firstCompDG = Kokkos::subview(dg, Kokkos::ALL(), std::make_pair(0, 1));
        Kokkos::deep_copy(firstCompDG, ma);
    }
}

/*!
 * @brief Copy data from an DGVector to a ModelArray through Kokkos views.
 *
 * @details Source and destination do not have to be located in the same memory space but strided
 * host -> device transfers can be extremely slow.
 *
 * @param ma The destination ModelArray view.
 * @param dg The source DGVector view.
 */
template <int N, typename... ArgsDG, typename... ArgsMA>
void kokkosDG2MA(const KokkosEigenView<ModelArray::DataType, ArgsMA...>& ma,
    const ConstKokkosEigenView<DGVector<N>, ArgsDG...>& dg)
{
    // N == 1 needs different treatment at compile-time because the corresponding Kokkos::View has
    // only one dimension.
    if constexpr (N == 1) {
        assert(ma.extent(1) == 1);
        // views need to have the same rank so we squeeze 1-sized component dimension
        const auto firstCompMA = Kokkos::subview(ma, Kokkos::ALL(), 0);
        Kokkos::deep_copy(firstCompMA, dg);
    } else if (N == ma.extent(1)) {
        Kokkos::deep_copy(ma, dg);
    } else {
        assert(ma.extent(1) == 1);
        // Assign only to the 0th component. Use a range to keep the dimension.
        const auto firstCompDG = Kokkos::subview(dg, Kokkos::ALL(), std::make_pair(0, 1));
        Kokkos::deep_copy(ma, firstCompDG);
    }
}

}

#endif // KOKKOSMODELARRAY_HPP