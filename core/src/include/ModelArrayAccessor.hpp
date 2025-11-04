/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "ModelArrayStore.hpp"

namespace Nextsim {

template <const TextTag& fieldName, bool isReadWrite = RO> class ModelArrayAccessor;

template <const TextTag& fieldName> class ModelArrayAccessor<fieldName, RO> {
public:
    ModelArrayAccessor(ModelArrayStore& store)
        : target(store.getRO(fieldName))
    {
    }

    const ModelArray& getHostRO()
    {
        if (target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(target.hostView(), target.deviceView());
            target.syncState = SyncState::SYNCED;
        }
        return target.modelArray;
    }

#ifdef USE_KOKKOS
    // returns a copy because target.deviceView has mutable data
    ConstDeviceView getDeviceRO()
    {
        DeviceView deviceView = target.deviceView();
        if (target.syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(deviceView, target.hostView());
            target.syncState = SyncState::SYNCED;
        }

        return deviceView;
    }
#endif

private:
    // lifetime and persistent address are enforced by ModelArrayStore
    ModelArrayStore::ExtModelArray& target;
};

template <const TextTag& fieldName> class ModelArrayAccessor<fieldName, RW> {
public:
    ModelArrayAccessor(ModelArrayStore& store)
        : target(store.getRW(fieldName))
    {
    }

    template <typename... Args>
    ModelArrayAccessor(ModelArrayStore& store, bool isReadWriteExternal, Args&&... args)
        : target(store._registerArray(fieldName, isReadWriteExternal, std::forward<Args>(args)...))
    {
    }

    const ModelArray& getHostRO()
    {
#ifdef USE_KOKKOS
        if (target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(target.hostView(), target.deviceView());
            target.syncState = SyncState::SYNCED;
        }
#endif
        return target.modelArray;
    }

    ModelArray& getHostRW()
    {
#ifdef USE_KOKKOS
        if (target.syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(target.hostView(), target.deviceView());
        }
        target.syncState = SyncState::HOST_CHANGED;
#endif

        return target.modelArray;
    }

#ifdef USE_KOKKOS
    // returns a copy because target.deviceView has mutable data
    ConstDeviceView getDeviceRO()
    {
        DeviceView deviceView = target.deviceView();
        if (target.syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(deviceView, target.hostView());
            target.syncState = SyncState::SYNCED;
        }

        return deviceView;
    }

    const DeviceView& getDeviceRW()
    {
        DeviceView deviceView = target.deviceView();
        if (target.syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(deviceView, target.hostView());
        target.syncState = SyncState::DEVICE_CHANGED;

        return deviceView;
    }
#endif

private:
    // lifetime and persistent address are enforced by ModelArrayStore
    ModelArrayStore::ExtModelArray& target;
};
}
