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

template <const TextTag& fieldName, bool isReadWrite = RO> class ModelArrayAccessor;

template <const TextTag& fieldName> class ModelArrayAccessor<fieldName, RO> {
public:
    ModelArrayAccessor(ModelArrayStore& store)
        : target(store.getRO(fieldName))
    {
    }

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

#ifdef USE_KOKKOS
    template <typename ExecSpace> const ModelArray& getHostRO(ExecSpace execSpace) const
    {
        if (target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(execSpace, target.hostView(), target.deviceView());
            target.syncState = SyncState::SYNCED;
        }
        return target.modelArray;
    }

    // returns a copy because target.deviceView has mutable data
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

    const ConstModelArrayAuto& getAutoRO() const
    {
#ifdef USE_KOKKOS
        return getDeviceRO();
#else
        return getHostRO();
#endif
    }
    template <typename ExecSpace> const ConstModelArrayAuto& getAutoRO(ExecSpace execSpace) const
    {
#ifdef USE_KOKKOS
        return getDeviceRO(execSpace);
#else
        return getHostRO();
#endif
    }

protected:
    // for the RW version
    ModelArrayAccessor(ModelArrayStore::ExtModelArray& _target)
        : target(_target)
    {
    }
    // lifetime and persistent address are enforced by ModelArrayStore
    ModelArrayStore::ExtModelArray& target;
};

template <const TextTag& fieldName>
class ModelArrayAccessor<fieldName, RW> : public ModelArrayAccessor<fieldName, RO> {
    using Base = ModelArrayAccessor<fieldName, RO>;

public:
    ModelArrayAccessor(ModelArrayStore& store)
        : Base(store.getRW(fieldName))
    {
    }

    template <typename... Args>
    ModelArrayAccessor(ModelArrayStore& store, bool isReadWriteExternal, Args&&... args)
        // using ModelArrayAccessor<fieldName, RO> directly instead of Base here leads to a compiler
        // error in ModelArrayAccessor_test.cpp
        : Base(store.registerArray(fieldName, isReadWriteExternal, std::forward<Args>(args)...))
    {
    }

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

#ifdef USE_KOKKOS
    template <typename ExecSpace> ModelArray& getHostRW(ExecSpace execSpace)
    {
        if (this->target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(execSpace, this->target.hostView(), this->target.deviceView());
        }
        this->target.syncState = SyncState::HOST_CHANGED;

        return this->target.modelArray;
    }

    DeviceModelArray& getDeviceRW()
    {
        if (this->target.syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(this->target.deviceView(), this->target.hostView());
        this->target.syncState = SyncState::DEVICE_CHANGED;

        return this->target.deviceModelArray();
    }

    template <typename ExecSpace> DeviceModelArray& getDeviceRW(ExecSpace execSpace)
    {
        if (this->target.syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(execSpace, this->target.deviceView(), this->target.hostView());
        this->target.syncState = SyncState::DEVICE_CHANGED;

        return this->target.deviceModelArray();
    }
#endif

    ModelArrayAuto& getAutoRW()
    {
#ifdef USE_KOKKOS
        return getDeviceRW();
#else
        return getHostRW();
#endif
    }
    
    template <typename ExecSpace> ModelArrayAuto& getAutoRW(ExecSpace execSpace)
    {
#ifdef USE_KOKKOS
        return getDeviceRW(execSpace);
#else
        return getHostRW();
#endif
    }
};
}

#endif /* MODELARRAYACCESSOR_HPP */