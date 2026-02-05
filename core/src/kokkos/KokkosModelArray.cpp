/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosModelArray.hpp"

namespace Nextsim {

ConstDeviceViewMA ConstDeviceModelArray::deviceView() const { return m_deviceView; }
ConstDeviceModelArray::operator ConstDeviceViewMA() const { return m_deviceView; }

const DeviceViewMA& DeviceModelArray::deviceView() { return m_deviceView; }
DeviceModelArray::operator const DeviceViewMA&() { return m_deviceView; }

void DeviceModelArray::assignData(const ConstDeviceModelArray& source)
{
    const ConstDeviceViewMA& srcView = source;
    if (m_deviceView.extent(0) != srcView.extent(0)) {
        const auto dstSubView = Kokkos::subview(m_deviceView, Kokkos::ALL(), std::make_pair(0, 1));
        const auto srcSubView = Kokkos::subview(srcView, Kokkos::ALL(), std::make_pair(0, 1));
        Kokkos::deep_copy(dstSubView, srcSubView);
    } else {
        Kokkos::deep_copy(m_deviceView, srcView);
    }
}

}