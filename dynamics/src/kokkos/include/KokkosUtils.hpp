/*!
 * @file KokkosUtils.hpp
 *
 * @date Feb 2, 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef KOKKOSUTILS_HPP
#define KOKKOSUTILS_HPP

#include <Eigen/Core>
#include <Kokkos_Bitset.hpp>
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

    template <typename Scalar> struct ToKokkosArrayDec<Scalar, -1, 1> {
        using Type = Scalar*;
    };

    template <typename Scalar> struct ToKokkosArrayDec<Scalar, 1, -1> {
        using Type = Scalar*;
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
    = Kokkos::View<typename Details::ToKokkosArrayDec<const typename EigenMat::Scalar,
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

/*!
 * @brief Creates compatible device and host views for a given Eigen matrix.
 */
template <typename EigenMat>
auto makeKokkosDualView(const std::string& name, EigenMat& mat, bool copy = false)
{
    return std::make_pair(makeKokkosHostView(mat), makeKokkosDeviceView(name, mat, copy));
}

template <typename T> using KokkosDeviceMapView = Kokkos::View<const T*>;

/*!
 * @brief Creates a const device view from an std::vector of simple data.
 *
 * @details The type T needs to be effectively trivially copyable, i.e. has no reference or pointer
 * members or non-default copy constructors. The caller has to ensure this, since enforcing the
 * requirements for T is impractical without C++ 20. The function works with compile-time-sized
 * Eigen matrices.
 *
 * @param name The name of the view.
 * @param buf The host side std::vector holding the data.
 * @param copy If true, the contents of mat are copied to the returned device buffer.
 */
template <typename T, typename Alloc>
auto makeKokkosDeviceViewMap(
    const std::string& name, const std::vector<T, Alloc>& buf, bool copy = false)
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

using DeviceBitset = Kokkos::Bitset<Kokkos::DefaultExecutionSpace>;
using ConstDeviceBitset = Kokkos::ConstBitset<Kokkos::DefaultExecutionSpace>;
/*!
 * @brief Creates a device bitset from an std::vector<bool>.
 *
 * @details The resulting Bitset is mutable but it can be assigned to a ConstBitset.
 *
 * @param buf The host side std::vector holding the data.
 */
DeviceBitset makeKokkosDeviceBitset(const std::vector<bool>& buf);

namespace Details {
    // Map Kokkos layout to Eigen layout options
    template <typename KokkosLayout> struct ToEigenLayout;

    template <> struct ToEigenLayout<Kokkos::LayoutLeft> {
        static constexpr int Options = Eigen::ColMajor;
    };
    template <> struct ToEigenLayout<Kokkos::LayoutRight> {
        static constexpr int Options = Eigen::RowMajor;
    };

    template <typename Scalar, int Rows, int Cols, int Options> struct ToMaybeConstMat {
        using Type = Eigen::Matrix<Scalar, Rows, Cols, Options>;
    };

    template <typename Scalar, int Rows, int Cols, int Options>
    struct ToMaybeConstMat<const Scalar, Rows, Cols, Options> {
        using Type = const Eigen::Matrix<Scalar, Rows, Cols, Options>;
    };

    // map kokkos array declaration to Eigen matrix
    template <class DataType, class Layout> struct ToEigenMatrix;

    template <typename Scalar, class Layout> struct ToEigenMatrix<Scalar**, Layout> {
        using Type = typename ToMaybeConstMat<Scalar, Eigen::Dynamic, Eigen::Dynamic,
            ToEigenLayout<Layout>::Options>::Type;
    };

    template <typename Scalar, int Cols, class Layout>
    struct ToEigenMatrix<Scalar* [Cols], Layout> {
        using Type = typename ToMaybeConstMat<Scalar, Eigen::Dynamic, Cols,
            (Cols == 1) ? Eigen::ColMajor : ToEigenLayout<Layout>::Options>::Type;
    };

    template <typename Scalar, int Rows, int Cols, class Layout>
    struct ToEigenMatrix<Scalar[Rows][Cols], Layout> {
        using Type = typename ToMaybeConstMat<Scalar, Rows, Cols,
            (Cols == 1) ? Eigen::ColMajor : ToEigenLayout<Layout>::Options>::Type;
    };

    template <typename Scalar, class Layout> struct ToEigenMatrix<Scalar*, Layout> {
        using Type = typename ToMaybeConstMat<Scalar, Eigen::Dynamic, 1, Eigen::ColMajor>::Type;
    };

    // map kokkos view spec to Eigen map
    template <class DataType, class Layout> struct ToEigenMap {
        using Type = Eigen::Map<typename ToEigenMatrix<DataType, Layout>::Type>;
    };

    /*    template <class DataType, class Layout> struct ToEigenMap<const
       std::remove_pointer_t<std::remove_all_extents_t<DataType>>, Layout> { using Type =
       Eigen::Map<const typename Details::ToEigenMatrix<DataType, Layout>::Type>;
        };*/
}

/*!
 * @brief Create an aquivalent Eigen map for a given Kokkos view.
 *
 * @details Supports 1D and 2D views. Constness of the underlying data is preserved.
 *
 * @param view The Kokkos view.
 */
template <class DataType, class... Properties>
KOKKOS_IMPL_FUNCTION auto makeEigenMap(const Kokkos::View<DataType, Properties...>& view)
{
    using View = Kokkos::View<DataType, Properties...>;
    static_assert(View::rank() <= 2, "Only 1D and 2D views can be converted to an Eigen map.");
    using MapType = typename Details::ToEigenMap<DataType, typename View::array_layout>::Type;

    if constexpr (View::rank() == 1) {
        return MapType(view.data(), view.extent(0));
    } else {
        return MapType(view.data(), view.extent(0), view.extent(1));
    }
}

} // namespace nextsim

#endif /* KOKKOSUTILS_HPP */