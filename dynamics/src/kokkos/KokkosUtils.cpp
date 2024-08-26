/*!
 * @file KokkosUtils.cpp
 *
 * @date Aug 26, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/KokkosUtils.hpp"

namespace Nextsim {

DeviceBitset makeKokkosDeviceBitset(const std::vector<bool>& buf)
{
    // unfortunately there is no more direct way to initialize a Kokkos::Bitset
    const unsigned nBits = buf.size();
    Kokkos::Bitset<Kokkos::HostSpace> bitsetHost(nBits);
    // fill with data
    bitsetHost.clear();
    for (unsigned i = 0; i < nBits; ++i) {
        if (buf[i])
            bitsetHost.set(i);
    }
    // host -> device
    Kokkos::Bitset<Kokkos::DefaultExecutionSpace> bitsetDevice(nBits);
    Kokkos::deep_copy(bitsetDevice, bitsetHost);

    return bitsetDevice;
}

} 