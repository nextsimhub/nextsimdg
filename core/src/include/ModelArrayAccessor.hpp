/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef MODELARRAYACCESSOR_HPP
#define MODELARRAYACCESSOR_HPP

#include "ModelArrayStore.hpp"

namespace Nextsim {

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

    template <typename ExecSpace> const ModelArray& getHostRO(ExecSpace execSpace) const
    {
#ifdef USE_KOKKOS
        if (target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(execSpace, target.hostView(), target.deviceView());
            target.syncState = SyncState::SYNCED;
        }
#endif
        return target.modelArray;
    }

#ifdef USE_KOKKOS
    // returns a copy because target.deviceView has mutable data
    ConstDeviceViewMA getDeviceRO() const
    {
        assert(target.modelArray.trueSize() > 0 && "ModelArray is allocated");

        DeviceViewMA deviceView = target.deviceView();
        if (target.syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(deviceView, target.hostView());
            target.syncState = SyncState::SYNCED;
        }

        return deviceView;
    }

    template <typename ExecSpace> ConstDeviceViewMA getDeviceRO(ExecSpace execSpace) const
    {
        assert(target.modelArray.trueSize() > 0 && "ModelArray is allocated");

        DeviceViewMA deviceView = target.deviceView();
        if (target.syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(execSpace, deviceView, target.hostView());
            target.syncState = SyncState::SYNCED;
        }

        return deviceView;
    }
#endif
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

    template <typename ExecSpace> ModelArray& getHostRW(ExecSpace execSpace)
    {
#ifdef USE_KOKKOS
        if (this->target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(execSpace, this->target.hostView(), this->target.deviceView());
        }
        this->target.syncState = SyncState::HOST_CHANGED;
#endif

        return this->target.modelArray;
    }

#ifdef USE_KOKKOS
    const DeviceViewMA& getDeviceRW()
    {
        const DeviceViewMA& deviceView = this->target.deviceView();
        if (this->target.syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(deviceView, this->target.hostView());
        this->target.syncState = SyncState::DEVICE_CHANGED;

        return deviceView;
    }

    template <typename ExecSpace> const DeviceViewMA& getDeviceRW(ExecSpace execSpace)
    {
        const DeviceViewMA& deviceView = this->target.deviceView();
        if (this->target.syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(execSpace, deviceView, this->target.hostView());
        this->target.syncState = SyncState::DEVICE_CHANGED;

        return deviceView;
    }
#endif

private:
    // lifetime and persistent address are enforced by ModelArrayStore
    //   ModelArrayStore::ExtModelArray& target;
};
}

#endif /* MODELARRAYACCESSOR_HPP */