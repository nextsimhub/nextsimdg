/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

// by also guarding for USE_KOKKOS this header can be safely included even
// when Kokkos is not enabled
#if !defined(KOKKOSMODELARRAY_HPP) && defined(USE_KOKKOS)
#define KOKKOSMODELARRAY_HPP

#include "../../include/ModelArray.hpp"
#include "KokkosUtils.hpp"
#include "dgVector.hpp"

namespace Nextsim {

// kokkos views compatible with ModelArray
using DeviceViewMA = KokkosDeviceView<ModelArray::DataType>;
using ConstDeviceViewMA = ConstKokkosDeviceView<ModelArray::DataType>;
using HostViewMA = KokkosHostView<ModelArray::DataType>;
using ConstHostViewMA = ConstKokkosHostView<ModelArray::DataType>;

template <int N>
void kokkosMA2DG(const ConstDeviceViewMA& ma, const KokkosDeviceView<DGVector<N>>& dg)
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

// template <int N> static ModelArray& kokkosDG2MA(const DGVector<N>& dg, ModelArray& ma)
// {
//     if (N == ma.components(0).size()) {
//         ma.setData(dg.data());
//     } else {
//         /* Assign the zero component as data. Since the setData function
//          * takes a pointer to continuous data, the data needs to be copied
//          * from the DGVector initially.
//          */
// ma.setData(dg.col(0));
// }
// return ma;
// }

}

#endif // KOKKOSMODELARRAY_HPP