/*!
 * @file KokkosUtils.hpp
 *
 * @date Feb 2, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSUTILS_HPP
#define KOKKOSUTILS_HPP

#include <Eigen/Core>
#include <Kokkos_Core.hpp>

#include <vector>

namespace Nextsim {

namespace Details {
    // Map Eigen matrix template parameters to the equivalent kokkos array type declaration.
    template <typename Scalar, int Rows, int Cols> struct ToKokkosArrayDec {
        using Type = Scalar[Rows][Cols];
    };

    template <typename Scalar> struct ToKokkosArrayDec<Scalar, -1, -1> {
        using Type = Scalar**;
    };

    template <typename Scalar, int Cols> struct ToKokkosArrayDec<Scalar, -1, Cols> {
        using Type = Scalar* [Cols];
    };

    // Map Eigen options to a kokkos layout.
    template <int Options> struct ToKokkosLayout {
        using Type = Kokkos::LayoutLeft;
    };
    template <> struct ToKokkosLayout<Eigen::RowMajor> {
        using Type = Kokkos::LayoutRight;
    };

}

// We can't specialize just for Eigen::Matrix because it needs to work with classes inheriting from
// Eigen::Matrix
template <typename EigenMat>
using KokkosDeviceView
    = Kokkos::View<typename Details::ToKokkosArrayDec<typename EigenMat::Scalar,
                       EigenMat::RowsAtCompileTime, EigenMat::ColsAtCompileTime>::Type,
        typename Details::ToKokkosLayout<EigenMat::Options>::Type>;
template <typename EigenMat>
using ConstKokkosDeviceView
    = Kokkos::View<const typename Details::ToKokkosArrayDec<typename EigenMat::Scalar,
                       EigenMat::RowsAtCompileTime, EigenMat::ColsAtCompileTime>::Type,
        typename Details::ToKokkosLayout<EigenMat::Options>::Type>;
template <typename EigenMat>
using KokkosHostView =
    typename Kokkos::View<typename Details::ToKokkosArrayDec<typename EigenMat::Scalar,
                              EigenMat::RowsAtCompileTime, EigenMat::ColsAtCompileTime>::Type,
        typename Details::ToKokkosLayout<EigenMat::Options>::Type, Kokkos::HostSpace,
        Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

/*!
 * @brief Creates a host view compatible with an Eigen matrix.
 *
 * @details The view points to the same memory as mat so it lifetime should be no longer than that
 * of mat.
 *
 * @param mat The Eigen matrix to use.
 */
template <typename EigenMat> auto makeKokkosHostView(EigenMat& mat)
{
    return KokkosHostView<EigenMat>(
        const_cast<typename EigenMat::Scalar*>(mat.data()), mat.rows(), mat.cols());
}

/*!
 * @brief Creates a device view compatible with an Eigen matrix.
 *
 * @details When building without device support a host view on mat is returned and no buffer is
 * allocated.
 *
 * @param mat The Eigen matrix to use.
 * @param copy If true, the contents of mat are copied to the returned device buffer.
 */
template <typename EigenMat>
auto makeKokkosDeviceView(const std::string& name, EigenMat& mat, bool copy = false)
{
    if constexpr (std::is_same_v<typename KokkosDeviceView<EigenMat>::memory_space,
                      Kokkos::HostSpace>) {
        return makeKokkosHostView(mat);
    } else {
        auto deviceView = KokkosDeviceView<EigenMat>(
            Kokkos::ViewAllocateWithoutInitializing(name), mat.rows(), mat.cols());

        if (copy) {
            auto hostView = makeKokkosHostView(mat);
            Kokkos::deep_copy(deviceView, hostView);
        }

        return deviceView;
    }
}

template <typename T> using KokkosDeviceMapView = Kokkos::View<const T*>;

template <typename T, typename Alloc>
auto makeKokkosDeviceViewMap(const std::string& name, const std::vector<T, Alloc>& buf, bool copy = false)
{
    using MapViewHost
        = Kokkos::View<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    if constexpr (std::is_same_v<typename KokkosDeviceMapView<T>::memory_space,
                      Kokkos::HostSpace>) {
        return MapViewHost(buf.data(), buf.size());
    } else {
        auto deviceView
            = Kokkos::View<T*>(Kokkos::ViewAllocateWithoutInitializing(name), buf.size());

        if (copy) {
            auto hostView = MapViewHost(buf.data(), buf.size());
            Kokkos::deep_copy(deviceView, hostView);
        }

        return KokkosDeviceMapView<T>(deviceView);
    }
}

}

#endif /* KOKKOSVPCGDYNAMICSKERNEL_HPP */