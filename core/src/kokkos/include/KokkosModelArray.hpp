/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

// by also guarding for USE_KOKKOS this header can be safely included even
// when Kokkos is not enabled
#if !defined(KOKKOSMODELARRAY_HPP) && defined(USE_KOKKOS)
#define KOKKOSMODELARRAY_HPP

#include "../../include/ModelArray.hpp"
#include "KokkosUtils.hpp"

namespace Nextsim {

// kokkos views compatible with ModelArray
using DeviceViewMA = KokkosDeviceView<ModelArray::DataType>;
using ConstDeviceViewMA = ConstKokkosDeviceView<ModelArray::DataType>;
using HostViewMA = KokkosHostView<ModelArray::DataType>;
using ConstHostViewMA = ConstKokkosHostView<ModelArray::DataType>;

}

#endif // KOKKOSMODELARRAY_HPP