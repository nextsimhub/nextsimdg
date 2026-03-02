/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef MODELARRAYACCESSOR_HPP
#define MODELARRAYACCESSOR_HPP

#include "ModelArrayStore.hpp"

namespace Nextsim {

// ModelArray type that is either on host or device, depending on the build options.
// See KokkosModelArray.hpp for why an additional type ConstModelArrayAuto is needed.
// Adding another const, i.e. "const ConstModelArrayAuto&", is possible but has no effect.
#ifdef USE_KOKKOS
using ModelArrayAuto = DeviceModelArray;
using ConstModelArrayAuto = ConstDeviceModelArray;
#else
using ModelArrayAuto = ModelArray;
// ModelArray has value semantics. The const applies to the data.
using ConstModelArrayAuto = const ModelArray;
#endif

// A handle that controls accesses through the underlying field
template <bool isReadWrite = RO> class ModelArrayAccessorBase;
// Variant with a compile-time name
template <const TextTag& fieldName, bool isReadWrite = RO> class ModelArrayAccessor;

// read-only implementation
template <> class ModelArrayAccessorBase<RO> {
public:
    /*!
     * @brief Construct a read-only accessor for a field in the given store.
     *
     * @details This does not create the actual field in the store. While the accessor may be
     * created beforehand, accesses to the field are only possible after it has been constructed
     * via the forwarding constructor of the read-write accessor.
     *
     * @param store The store in which the field is registered.
     * @param fieldName Identifier of the field in the store.
     */
    ModelArrayAccessorBase(ModelArrayStore& store, const std::string& fieldName)
        : target(store.getRO(fieldName))
    {
    }

    /// @brief Get read-only access to the field on the CPU.
    const ModelArray& getHostRO() const
    {
#ifdef USE_KOKKOS
        if (target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(target.hostView(), target.deviceView());
            target.syncState = SyncState::SYNCED;
        }
#endif
        return target.modelArray;
    }

    /*!
     * @brief Get read-only access to the field on the CPU (asynchronous).
     *
     * @details With Kokkos enabled, the provided execution space is used for host-device data
     * transfers. This results in an asynchronous copy that is not necessarily finished by the time
     * this function returns.
     *
     * @param execSpace The Kokkos execution space or a dummy object.
     */
#ifdef USE_KOKKOS
    template <typename ExecSpace> const ModelArray& getHostRO(ExecSpace execSpace) const
    {
        if (target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(execSpace, target.hostView(), target.deviceView());
            target.syncState = SyncState::SYNCED;
        }
        return target.modelArray;
    }

    /// @brief Get read-only access to the field on the GPU.
    const ConstDeviceModelArray& getDeviceRO() const
    {
        assert(target.modelArray.trueSize() > 0 && "ModelArray is allocated");

        DeviceViewMA deviceView = target.deviceView();
        if (target.syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(target.deviceView(), target.hostView());
            target.syncState = SyncState::SYNCED;
        }

        return target.deviceModelArray();
    }

    /*!
     * @brief Get read-only access to the field on the GPU (asynchronous).
     *
     * @details With Kokkos enabled, the provided execution space is used for host-device data
     * transfers. This results in an asynchronous copy that is not necessarily finished by the time
     * this function returns.
     *
     * @param execSpace The Kokkos execution space or a dummy object.
     */
    template <typename ExecSpace>
    const ConstDeviceModelArray& getDeviceRO(ExecSpace execSpace) const
    {
        assert(target.modelArray.trueSize() > 0 && "ModelArray is allocated");

        if (target.syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(execSpace, target.deviceView(), target.hostView());
            target.syncState = SyncState::SYNCED;
        }

        return target.deviceModelArray();
    }
#endif

    /// @brief Get read-only access to the field on the GPU (Kokkos enabled) or the CPU.
    const ConstModelArrayAuto& getAutoRO() const
    {
#ifdef USE_KOKKOS
        return getDeviceRO();
#else
        return getHostRO();
#endif
    }

    /// @brief Get read-only access to the field on the GPU (Kokkos enabled) or the CPU
    /// (asynchronous).
    template <typename ExecSpace> const ConstModelArrayAuto& getAutoRO(ExecSpace execSpace) const
    {
#ifdef USE_KOKKOS
        return getDeviceRO(execSpace);
#else
        return getHostRO();
#endif
    }

    /*!
     * @brief Collects read-only accessors for every field in the given store.
     *
     * @param store The model array store to read from.
     * @return A map of field names to accessors.
     */
    static std::unordered_map<std::string, ModelArrayAccessorBase<RO>> getAll(
        const ModelArrayStore& store)
    {
        std::unordered_map<std::string, ModelArrayAccessorBase<RO>> dataMap;

        for (const auto& [name, extArrFlagged] : store.store) {
            // skip fields that don't actually exist
            if (!extArrFlagged.isRegistered) {
                continue;
            }
            // Internally ModelArrayAccessor always holds a mutable reference but the RO variant
            // only exposes the data as const so this is safe.
            dataMap.emplace(name,
                ModelArrayAccessorBase<RO>(
                    const_cast<ModelArrayStore::ExtModelArray&>(extArrFlagged.extModelArray)));
        }

        return dataMap;
    }

protected:
    // for the RW version
    ModelArrayAccessorBase(ModelArrayStore::ExtModelArray& _target)
        : target(_target)
    {
    }
    // lifetime and persistent address are enforced by ModelArrayStore
    ModelArrayStore::ExtModelArray& target;
};

// read-only implementation with added compile-time field name
template <const TextTag& fieldName>
class ModelArrayAccessor<fieldName, RO> : public ModelArrayAccessorBase<RO> {
public:
    ModelArrayAccessor(ModelArrayStore& store)
        : ModelArrayAccessorBase<RO>(store, fieldName)
    {
    }
};

// read-write implementation
template <> class ModelArrayAccessorBase<RW> : public ModelArrayAccessorBase<RO> {
    using Base = ModelArrayAccessorBase<RO>;

public:
    /*!
     * @brief Construct a read-write accessor for a field in the given store.
     *
     * @details This does not create the actual field in the store. While the accessor may be
     * created beforehand, accesses to the field are only possible after it has been constructed
     * via the forwarding constructor of the read-write accessor.
     *
     * @param store The store in which the field is registered.
     * @param fieldName Identifier of the field in the store.
     */
    ModelArrayAccessorBase(ModelArrayStore& store, const std::string& fieldName)
        : Base(store.getRW(fieldName))
    {
    }

    /*!
     * @brief Construct a field and a read-write accessor for it in the given store.
     *
     * @details For every field in a given store this constructor has to be used at least once to
     * actually allocate and construct the field. Repeated calls are allowed if isReadWriteExternal
     * is the same and will cause the field to be recreated inplace.
     *
     * @param store The store in which the field is registered.
     * @param fieldName Identifier of the field in the store.
     * @param isReadWriteExternal Determines whether non-constructing read-write accessors can be
     * created for this field.
     * @param args Arguments that are forwarded to the constructor of ModelArray.
     */
    template <typename... Args>
    ModelArrayAccessorBase(ModelArrayStore& store, const std::string& fieldName,
        bool isReadWriteExternal, Args&&... args)
        // using ModelArrayAccessor<fieldName, RO> directly instead of Base here leads to a compiler
        // error in ModelArrayAccessor_test.cpp
        : Base(store.registerArray(fieldName, isReadWriteExternal, std::forward<Args>(args)...))
    {
    }

    /// @brief Get read-only access to the field on the GPU (Kokkos enabled) or the CPU.
    ModelArray& getHostRW()
    {
#ifdef USE_KOKKOS
        if (this->target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(this->target.hostView(), this->target.deviceView());
        }
        this->target.syncState = SyncState::HOST_CHANGED;
#endif

        return this->target.modelArray;
    }

    /// @brief Get read-write access to the field on the CPU (asynchronous).
#ifdef USE_KOKKOS
    template <typename ExecSpace> ModelArray& getHostRW(ExecSpace execSpace)
    {
        if (this->target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(execSpace, this->target.hostView(), this->target.deviceView());
        }
        this->target.syncState = SyncState::HOST_CHANGED;

        return this->target.modelArray;
    }

    /// @brief Get read-only access to the field on the GPU.
    DeviceModelArray& getDeviceRW()
    {
        if (this->target.syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(this->target.deviceView(), this->target.hostView());
        this->target.syncState = SyncState::DEVICE_CHANGED;

        return this->target.deviceModelArray();
    }

    /// @brief Get read-write access to the field on the GPU (asynchronous).
    template <typename ExecSpace> DeviceModelArray& getDeviceRW(ExecSpace execSpace)
    {
        if (this->target.syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(execSpace, this->target.deviceView(), this->target.hostView());
        this->target.syncState = SyncState::DEVICE_CHANGED;

        return this->target.deviceModelArray();
    }
#endif

    /// @brief Get read-write access to the field on the GPU (Kokkos enabled) or the CPU.
    ModelArrayAuto& getAutoRW()
    {
#ifdef USE_KOKKOS
        return getDeviceRW();
#else
        return getHostRW();
#endif
    }

    /// @brief Get read-only access to the field on the GPU (Kokkos enabled) or the CPU
    /// (asynchronous).
    template <typename ExecSpace> ModelArrayAuto& getAutoRW(ExecSpace execSpace)
    {
#ifdef USE_KOKKOS
        return getDeviceRW(execSpace);
#else
        return getHostRW();
#endif
    }
};

// read-write implementation with added compile-time field name
template <const TextTag& fieldName>
class ModelArrayAccessor<fieldName, RW> : public ModelArrayAccessorBase<RW> {
    using Base = ModelArrayAccessorBase<RW>;

public:
    ModelArrayAccessor(ModelArrayStore& store)
        : Base(store, fieldName)
    {
    }

    template <typename... Args>
    ModelArrayAccessor(ModelArrayStore& store, bool isReadWriteExternal, Args&&... args)
        : Base(store, fieldName, isReadWriteExternal, std::forward<Args>(args)...)
    {
    }
};

}

#endif /* MODELARRAYACCESSOR_HPP */