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

/* Wrapper for Kokkos views with semantics closer to ModelArray.
 * The fundamental difference is that ModelArray has value semantics, while a Kokkos view has
 * reference semantics. When accessing either through the ModelArrayStore the provided interface can
 * be similar since whe always work on a reference. However, constness has to be part of the type
 * here to ensure proper capture-by-value in lambdas for kernel execution. ConstDeviceModelArray is
 * essentially just a "const DeviceModelArray".
 */
class ConstDeviceModelArray {
public:
    KOKKOS_IMPL_FUNCTION double operator[](DeviceIndex i) const { return m_deviceView(i, 0); }

    // provides a copy because the underlying view is not const
    ConstDeviceViewMA deviceView() const;
    operator ConstDeviceViewMA() const;

protected:
    ConstDeviceModelArray() = default;

    // not const so that it can be reused in DeviceModelArray
    DeviceViewMA m_deviceView;

    friend class ModelArrayStore;
};
class DeviceModelArray : public ConstDeviceModelArray {
public:
    // Functions callable within a kernel should be const, even when they mutate the underlying
    // data. When captured by copy within a kernel lambda, any DeviceModelArray will be constant by
    // default. To prevent this, every lambda would need to be marked as mutable.
    KOKKOS_IMPL_FUNCTION double& operator[](DeviceIndex i) const { return m_deviceView(i, 0); }

    const DeviceViewMA& deviceView();

    operator const DeviceViewMA&();

    void assignData(const ConstDeviceModelArray& source);

private:
    DeviceModelArray() = default;

    friend class ModelArrayStore;
};

}

#endif // KOKKOSMODELARRAY_HPP