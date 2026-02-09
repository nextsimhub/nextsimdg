/*!
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

// by also guarding for USE_KOKKOS this header can be safely included even
// when Kokkos is not enabled
#if !defined(KOKKOSUTILS_HPP) && defined(USE_KOKKOS)
#define KOKKOSUTILS_HPP

#include <Eigen/Core>
#include <Kokkos_Bitset.hpp>
#include <Kokkos_Core.hpp>

#include <vector>

namespace Nextsim {

using DeviceIndex = EIGEN_DEFAULT_DENSE_INDEX_TYPE;

template <typename ExecutionSpace> constexpr inline bool IS_GPU_EXEC_SPACE = false;

#ifdef KOKKOS_ENABLE_CUDA
template <> constexpr inline bool IS_GPU_EXEC_SPACE<Kokkos::Cuda> = true;
#endif

#ifdef KOKKOS_ENABLE_HIP
template <> constexpr inline bool IS_GPU_EXEC_SPACE<Kokkos::HIP> = true;
#endif

// Land checks currently only improve performance on GPU. On CPU they lead to a slowdown so checks
// should only be introduced when necessary. Currently there are no optional land checks.
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
// We can't directly use if constexpr(IS_GPU_EXEC_SPACE) in kernels because lambda capture for
// values only appearing in constexpr branches is not supported in CUDA.
static_assert(IS_GPU_EXEC_SPACE<Kokkos::DefaultExecutionSpace>, "expected gpu execution space");
#define ENABLE_OPTIONAL_LAND_CHECKS
#endif

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

// Kokkos view types that are compatible (type, compile-time size, storage order, constness) with
// Eigen matrices. We can't specialize just for Eigen::Matrix because it needs to work with classes
// inheriting from Eigen::Matrix like DGVector
template <typename EigenMat, typename... Args>
using KokkosEigenView
    = Kokkos::View<typename Details::ToKokkosArrayDec<typename EigenMat::Scalar,
                       EigenMat::RowsAtCompileTime, EigenMat::ColsAtCompileTime>::Type,
        typename Details::ToKokkosLayout<EigenMat::Options>::Type, Args...>;
template <typename EigenMat, typename... Args>
using ConstKokkosEigenView
    = Kokkos::View<typename Details::ToKokkosArrayDec<const typename EigenMat::Scalar,
                       EigenMat::RowsAtCompileTime, EigenMat::ColsAtCompileTime>::Type,
        typename Details::ToKokkosLayout<EigenMat::Options>::Type>;

template <typename EigenMat> using KokkosDeviceView = KokkosEigenView<EigenMat>;
template <typename EigenMat> using ConstKokkosDeviceView = ConstKokkosEigenView<EigenMat>;
template <typename EigenMat>
using KokkosHostView
    = KokkosEigenView<EigenMat, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
template <typename EigenMat>
using ConstKokkosHostView
    = ConstKokkosEigenView<EigenMat, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

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
    // const_cast is necessary because Eigen only gives const access to the underlying pointer
    if constexpr (EigenMat::RowsAtCompileTime == 1 || EigenMat::ColsAtCompileTime == 1) {
        return KokkosHostView<EigenMat>(
            const_cast<typename EigenMat::Scalar*>(mat.data()), mat.rows() * mat.cols());
    }
    return KokkosHostView<EigenMat>(
        const_cast<typename EigenMat::Scalar*>(mat.data()), mat.rows(), mat.cols());
}
// const overload
template <typename EigenMat> auto makeKokkosHostView(const EigenMat& mat)
{
    if constexpr (EigenMat::RowsAtCompileTime == 1 || EigenMat::ColsAtCompileTime == 1) {
        return ConstKokkosHostView<EigenMat>(mat.data(), mat.rows() * mat.cols());
    }
    return ConstKokkosHostView<EigenMat>(mat.data(), mat.rows(), mat.cols());
}

/// Options for the creation of views based on existing data.
enum struct MakeViewOptions {
    ALWAYS_COPY, ///< The resulting view always owns its buffer and it is initialized with the given
                 ///< data. Use when the created view needs to outlive the source data.
    DEVICE_COPY, ///< Copy the data only when it is not in the right memory space. Host views
                 ///< created with this option will just point to the existing data.
    NO_COPY ///< Perform no initialization.
};

/*!
 * @brief Creates a device view compatible with an Eigen matrix.
 *
 * @details When building without device support a host view on mat is returned and a buffer is
 * only allocated with MakeViewOptions::ALWAYS_COPY.
 *
 * @param mat The Eigen matrix to use.
 * @param opts See MakeViewOptions.
 */
template <typename EigenMat>
auto makeKokkosDeviceView(
    const std::string& name, EigenMat& mat, MakeViewOptions opts = MakeViewOptions::NO_COPY)
{
    if constexpr (std::is_same_v<typename KokkosDeviceView<EigenMat>::memory_space,
                      Kokkos::HostSpace>) {
        return makeKokkosHostView(mat);
    } else {
        auto deviceView = [&]() {
            // 1D matrix. Using a two arg constructor works in both cases but kokkos
            // complains when debugging is enabled.
            if constexpr (EigenMat::RowsAtCompileTime == 1 || EigenMat::ColsAtCompileTime == 1) {
                return KokkosDeviceView<EigenMat>(
                    Kokkos::ViewAllocateWithoutInitializing(name), mat.rows() * mat.cols());
            } else {
                return KokkosDeviceView<EigenMat>(
                    Kokkos::ViewAllocateWithoutInitializing(name), mat.rows(), mat.cols());
            }
        }();

        if (opts != MakeViewOptions::NO_COPY) {
            auto hostView = makeKokkosHostView(mat);
            Kokkos::deep_copy(deviceView, hostView);
        }

        return deviceView;
    }
}

/*!
 * @brief Creates compatible device and host views for a given Eigen matrix.
 */
template <typename EigenMat> auto makeKokkosDualView(const std::string& name, EigenMat& mat)
{
    return std::make_pair(
        makeKokkosHostView(mat), makeKokkosDeviceView(name, mat, MakeViewOptions::NO_COPY));
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
 * @param opts See MakeViewOptions.
 */
template <typename T, typename Alloc>
auto makeKokkosDeviceViewMap(const std::string& name, const std::vector<T, Alloc>& buf,
    MakeViewOptions opts = MakeViewOptions::NO_COPY)
{
    using MapViewHost
        = Kokkos::View<const T*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

    if constexpr (std::is_same_v<typename KokkosDeviceMapView<T>::memory_space,
                      Kokkos::HostSpace>) {
        if (opts == MakeViewOptions::ALWAYS_COPY) {
            // need a mutable view to copy the data
            auto view = Kokkos::View<T*>(Kokkos::ViewAllocateWithoutInitializing(name), buf.size());
            auto unmanagedView = MapViewHost(buf.data(), buf.size());

            Kokkos::deep_copy(view, unmanagedView);
            // create new view to make it immutable
            return KokkosDeviceMapView<T>(view);
        }
        return KokkosDeviceMapView<T>(MapViewHost(buf.data(), buf.size()));
    } else {
        auto deviceView
            = Kokkos::View<T*>(Kokkos::ViewAllocateWithoutInitializing(name), buf.size());

        if (opts != MakeViewOptions::NO_COPY) {
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
 * @brief Create an equivalent Eigen map for a given Kokkos view.
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

template <typename T> constexpr KOKKOS_IMPL_FUNCTION T sqr(T x) { return x * x; }

} // namespace nextsim

#endif /* KOKKOSUTILS_HPP */