/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef MASTORE_HPP
#define MASTORE_HPP

#include "../kokkos/include/KokkosModelArray.hpp"
#include "ModelArray.hpp"
#include "ModelArrayRef.hpp" // for RW,RO globals

#include <string>
#include <unordered_map>

struct TextTag;

namespace Nextsim {

#ifdef USE_KOKKOS

// Wrapper for Kokkos views with semantics closer to ModelArray
class ConstDeviceModelArray {
public:
    ConstDeviceModelArray(ConstDeviceViewMA deviceView)
        : m_deviceView(std::move(deviceView))
    {
    }

    KOKKOS_IMPL_FUNCTION double operator[](DeviceIndex i) const { return m_deviceView(i, 0); }

    operator ConstDeviceViewMA() const;

private:
    ConstDeviceViewMA m_deviceView;

    friend class ModelArrayStore;
};
class DeviceModelArray {
public:
    // DeviceModelArray(DeviceViewMA deviceView) : m_deviceView(std::move(m_deviceView)) {}

    KOKKOS_IMPL_FUNCTION double& operator[](DeviceIndex i) const { return m_deviceView(i, 0); }

    operator const DeviceViewMA&() const;

private:
    DeviceModelArray() = default;
    DeviceViewMA m_deviceView;

    friend class ModelArrayStore;
};

enum struct SyncState { SYNCED, HOST_CHANGED, DEVICE_CHANGED };
#endif

class ModelArrayStore {
public:
    std::unordered_map<std::string, const ModelArray*> getAllData() const;

    // Checks that all known fields have been properly registered with an Accessor that provided the
    // constructor arguments.
    // @return True iff all fields are registered.
    bool checkAllRegistered() const;

private:
    struct ExtModelArray {
        std::string name;
        ModelArray modelArray;

#ifdef USE_KOKKOS
        SyncState syncState = SyncState::SYNCED;
        // not a simple data member because the host buffer owned by ModelArray can be overwritten
        HostViewMA hostView();
        // handles lazy initialization for the device buffer
        const DeviceViewMA& deviceView();
        DeviceModelArray& deviceModelArray();

    private:
        DeviceModelArray m_deviceModelArray;
#endif
    };

    struct ExtModelArrayFlagged {
        ExtModelArrayFlagged(const std::string& name, bool _isReadWrite, bool _isRegistered);

        bool isReadWrite;
        bool isRegistered;
        ExtModelArray extModelArray;
    };

    template <typename... Args>
    ExtModelArray& registerArray(const std::string& field, bool isReadWrite, Args&&... args)
    {
        auto it = store.find(field);
        if (it != store.end()) {
            ExtModelArrayFlagged& extArrFlagged = it->second;
            if (extArrFlagged.isReadWrite == RW && isReadWrite == RO) {
                throw std::logic_error("Registering ModelArray \"" + field
                    + "\" as read-only but at least one accessor requested read-write.");
            }

            if (extArrFlagged.isRegistered) {
                // Double registration unfortunately has to be allowed because two instances
                // of each ModelComponent are created which try to register the same fields.
                // This loophole can be used to get a RW access for a field that should be RO.
                // throw std::logic_error(
                //    "Registering ModelArray \"" + field + "\" but it was already registered.");
                if (extArrFlagged.isReadWrite != isReadWrite) {
                    throw std::logic_error("Registering ModelArray \"" + field
                        + "\" but it was already registered with different access restrictions.");
                }
            } else {
                extArrFlagged.isReadWrite = isReadWrite;
                extArrFlagged.isRegistered = true;
            }
        } else {
            // Regular emplace would be fine here since we know that it does not exist but
            // try_emplace has a more ergonomic signature.
            it = store.try_emplace(field, field, isReadWrite, true).first;
        }

        ExtModelArray& extArr = it->second.extModelArray;

        // destroy and construct inplace to keep pointers valid
        extArr.modelArray.~ModelArray();
        new (&extArr.modelArray) ModelArray(std::forward<Args>(args)...);

        return extArr;
    }

    ExtModelArray& getRW(const std::string& field);
    ExtModelArray& getRO(const std::string& field);

    std::unordered_map<std::string, ExtModelArrayFlagged> store;

    template <const TextTag& fieldName, bool isReadWrite> friend class ModelArrayAccessor;
};

}

#endif /* MARSTORE_HPP */
