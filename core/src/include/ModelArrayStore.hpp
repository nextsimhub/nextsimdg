/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef MASTORE_HPP
#define MASTORE_HPP

#include "../../../dynamics/src/kokkos/include/KokkosUtils.hpp"
#include "ModelArray.hpp"
#include "ModelArrayRef.hpp" // for RW,RO globals
#include "ModelArrayReferenceStore.hpp" // only for testing

#include <string>
#include <unordered_map>

struct TextTag;

namespace Nextsim {

// kokkos views compatible with ModelArray
using DeviceView = KokkosDeviceView<ModelArray::DataType>;
using ConstDeviceView = ConstKokkosDeviceView<ModelArray::DataType>;
using HostView = KokkosHostView<ModelArray::DataType>;
using ConstHostView = ConstKokkosHostView<ModelArray::DataType>;

// Inherit from ModelArrayReferenceStore only for demonstration purposes so we can replace
// ModelArrayReferenceStore with ModelArrayStore without changing everything.
class ModelArrayStore : public ModelArrayReferenceStore {
public:
    // Resizes all arrays and allocates device buffers.
    // Needs to be done after the ModelArray::m_sz was properly initialized and should be done only
    // once after all arrays have been registered.
    void resize_arrays();

private:
    struct ExtModelArray {
        ModelArray modelArray;
        DeviceView deviceView;

        // not a simple data member because the host buffer owned by ModelArray can be overwritten
        HostView hostView();
        ConstHostView hostView() const;
        // HostView hostView;
    };

    struct ExtModelArrayFlagged {
        ExtModelArrayFlagged(bool _isReadWrite, bool _isRegistered)
            : isReadWrite(_isReadWrite)
            , isRegistered(_isRegistered)
        {
        }

        bool isReadWrite;
        bool isRegistered;
        ExtModelArray extModelArray;
    };

    // underscore in name is just to prevent name conflict with ModelArrayRefStore for now
    template <typename... Args>
    ExtModelArray& _registerArray(const std::string& field, bool isReadWrite, Args&&... args)
    {
        if (was_resized) {
            throw std::logic_error(
                "Registering ModelArray \"" + field + "\" after resize_arrays() was called.");
        }

        auto it = store.find(field);
        if (it != store.end()) {
            ExtModelArrayFlagged& extArrFlagged = it->second;
            if (extArrFlagged.isReadWrite == RW && isReadWrite == RO) {
                throw std::logic_error("Registering ModelArray \"" + field
                    + "\" as read-only but at least one accessor requested read-write.");
            }
            // Double registration where previously allowed but it could cause unexpected behaviour
            // when the ModelArray params do not agree and should not be a need for it.
            if (extArrFlagged.isRegistered) {
                throw std::logic_error(
                    "Registering ModelArray \"" + field + "\" but it was already registered.");
            }
            extArrFlagged.isReadWrite = isReadWrite;
            extArrFlagged.isRegistered = true;
        } else {
            // Regular emplace would be fine here since we know that it does not exist but
            // try_emplace has a more ergonomic signature.
            it = store.try_emplace(field, isReadWrite, true).first;
        }

        ExtModelArray& extArr = it->second.extModelArray;

        // destroy and construct inplace to keep pointers valid
        extArr.modelArray.~ModelArray();
        new (&extArr.modelArray) ModelArray(std::forward<Args>(args)...);

        return extArr;
    }

    ExtModelArray& getRW(const std::string& field);
    const ExtModelArray& getRO(const std::string& field);

    std::unordered_map<std::string, ExtModelArrayFlagged> store;
    bool was_resized = false;

    template <const TextTag& fieldName, bool isReadWrite> friend class ModelArrayAccessor;
};

}

#endif /* MARSTORE_HPP */
