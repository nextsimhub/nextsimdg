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

namespace Nextsim {

#ifdef USE_KOKKOS
enum struct SyncState { SYNCED, HOST_CHANGED, DEVICE_CHANGED };
#endif

class ModelArrayStore {
public:
    /*!
     * @brief Checks that all known fields have been properly registered with an Accessor that
     * provided the constructor arguments.
     *
     * @details The purpose of this function is to detect a faulty initialization early. The
     * results can be ignored and may not indicate that the program does not work. However, a field
     * that was not properly registered is just default constructed and may have an unexpected type
     * or size that can cause subtle bugs later on. If a field is not used, it should not be
     * registered at all.
     *
     * @return True iff all fields are properly registered.
     */
    bool checkAllRegistered() const;

private:
    // extended ModelArray with compatible KokkosView and name
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

    // ModelArray with additional internal flags to track the registration state
    struct ExtModelArrayFlagged {
        ExtModelArrayFlagged(const std::string& name, bool _isReadWrite, bool _isRegistered);

        bool isReadWrite;
        bool isRegistered; // extModelArray was created through a call of registerArray
        ExtModelArray extModelArray;
    };

    // (re)create a ModelArray in the store
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

    // access to existing fields
    ExtModelArray& getRW(const std::string& field);
    ExtModelArray& getRO(const std::string& field);

    std::unordered_map<std::string, ExtModelArrayFlagged> store;

    // give access to the getters and registerArray
    template <bool isReadWrite> friend class ModelArrayAccessorBase;
};

}

#endif /* MARSTORE_HPP */
