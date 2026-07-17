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
using HostViewMA = KokkosHostViewUnmanaged<ModelArray::DataType>;
using ConstHostViewMA = ConstKokkosHostViewUnmanaged<ModelArray::DataType>;

/*!
 * @brief Wrapper for const Kokkos views with semantics closer to ModelArray for use in auto
 * contexts.
 *
 * @details The fundamental difference is that ModelArray has value semantics, while a Kokkos view
 * has reference semantics. When accessing either through the ModelArrayStore the provided interface
 * can be similar since whe always work on a reference. However, constness has to be part of the
 * type here to ensure proper capture-by-value in lambdas for kernel execution.
 * ConstDeviceModelArray is essentially just a "const DeviceModelArray".
 */
class ConstDeviceModelArray {
public:
    KOKKOS_IMPL_FUNCTION FloatType operator[](DeviceIndex i) const { return m_deviceView(i, 0); }

    // provides a copy because the underlying view is not const
    ConstDeviceViewMA deviceView() const;
    operator ConstDeviceViewMA() const;

protected:
    ConstDeviceModelArray() = default;

    // not const so that it can be reused in DeviceModelArray
    DeviceViewMA m_deviceView;

    friend class ModelArrayStore;
};

/// @brief Wrapper for mutable Kokkos views with semantics closer to ModelArray.
class DeviceModelArray : public ConstDeviceModelArray {
public:
    /*!
     * @brief Accessor for the 0-component of the field, similar to ModelArray::operator[].
     *
     * @details Functions callable within a kernel should be const, even when they mutate the
     * underlying data. When captured by copy within a kernel lambda, any DeviceModelArray will be
     * constant by default. To prevent this, every lambda would need to be marked as mutable.
     */
    //
    KOKKOS_IMPL_FUNCTION FloatType& operator[](DeviceIndex i) const { return m_deviceView(i, 0); }

    const DeviceViewMA& deviceView();
    operator const DeviceViewMA&();

    /// @brief Sets data from another DeviceModelArray, respecting the type of this, similar to
    /// ModelArray::assignData.
    void assignData(const ConstDeviceModelArray& source);

private:
    DeviceModelArray() = default;

    friend class ModelArrayStore;
};

}

#endif // KOKKOSMODELARRAY_HPP
