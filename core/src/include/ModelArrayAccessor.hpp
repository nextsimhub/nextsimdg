/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "ModelArrayStore.hpp"

namespace Nextsim {

enum struct SyncState { HOST_CHANGED, DEVICE_CHANGED, SYNCED };

template <const TextTag& fieldName, bool isReadWrite> class ModelArrayAccessor;

template <const TextTag& fieldName> class ModelArrayAccessor<fieldName, RO> {
public:
    ModelArrayAccessor(ModelArrayStore& store)
        : target(store.getRO(fieldName))
    {
    }

    const ModelArray& getHostRO()
    {
        if (syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(target.hostView(), target.deviceView);
            syncState = SyncState::SYNCED;
        }
        return target.modelArray;
    }

    // returns a copy because target.deviceView has mutable data
    ConstDeviceView getDeviceRO()
    {
        if (syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(target.deviceView, target.hostView());
            syncState = SyncState::SYNCED;
        }

        return target.deviceView;
    }

private:
    // mirrors DeviceView for data transfers but is not part of the interface
    using HostView = KokkosHostView<ModelArray::DataType>;

    SyncState syncState;

    // lifetime and persistent address are enforced by ModelArrayStore
    const ModelArrayStore::ExtModelArray& target;
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
    {}

    const ModelArray& getHostRO()
    {
        if (syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(target.hostView(), target.deviceView);
            syncState = SyncState::SYNCED;
        }
        return target.modelArray;
    }

    ModelArray& getHostRW()
    {
        if (syncState == SyncState::DEVICE_CHANGED) {
            Kokkos::deep_copy(target.hostView(), target.deviceView);
        }
        syncState = SyncState::HOST_CHANGED;

        return target.modelArray;
    }

    // returns a copy because target.deviceView has mutable data
    ConstDeviceView getDeviceRO()
    {
        if (syncState == SyncState::HOST_CHANGED) {
            Kokkos::deep_copy(target.deviceView, target.hostView());
            syncState = SyncState::SYNCED;
        }

        return target.deviceView;
    }

    const DeviceView& getDeviceRW()
    {
        if (syncState == SyncState::HOST_CHANGED)
            Kokkos::deep_copy(target.deviceView, target.hostView());
        syncState = SyncState::DEVICE_CHANGED;

        return target.deviceView;
    }

private:
    SyncState syncState;

    // lifetime and persistent address are enforced by ModelArrayStore
    ModelArrayStore::ExtModelArray& target;
};
}
