/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef MASTORE_HPP
#define MASTORE_HPP

#include "../../../dynamics/src/kokkos/include/KokkosUtils.hpp"
#include "ModelArray.hpp"
#include "ModelArrayReferenceStore.hpp" // only for testing

#include <string>
#include <unordered_map>

struct TextTag;

namespace Nextsim {

// kokkos views compatible with ModelArray
using DeviceView = KokkosDeviceView<ModelArray::DataType>;
using ConstDeviceView = ConstKokkosDeviceView<ModelArray::DataType>;
using HostView = KokkosHostView<ModelArray::DataType>;

// Inherit from ModelArrayReferenceStore only for demonstration purposes so we can replace
// ModelArrayReferenceStore with ModelArrayStore without changing everything.
class ModelArrayStore : public ModelArrayReferenceStore {

    struct ExtModelArray {
        ModelArray modelArray;
        DeviceView deviceView;
        HostView hostView;
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

        // initialize the views which need a pointer to the actual underlying buffer
        std::tie(extArr.hostView, extArr.deviceView)
            = makeKokkosDualView(field, extArr.modelArray.getDataRef());

        return extArr;
    }

    ExtModelArray& getRW(const std::string& field)
    {
        if (auto it = store.find(field); it != store.end()) {
            ExtModelArrayFlagged& extArrFlagged = it->second;
            if (!extArrFlagged.isReadWrite) {
                if (extArrFlagged.isRegistered) {
                    throw std::logic_error("Trying to access the read-only ModelArray \"" + field
                        + "\" as read-write.");
                } else {
                    // promote to read-write because we don't now the true restriction yet
                    extArrFlagged.isReadWrite = true;
                }
            }
            return extArrFlagged.extModelArray;
        }

        // Regular emplace would be fine here since we know that it does not exist but try_emplace
        // has a more ergonomic signature.
        return store.try_emplace(field, RW, false).first->second.extModelArray;
    }

    const ExtModelArray& getRO(const std::string& field)
    {
        return store.try_emplace(field, RO, false).first->second.extModelArray;
    }

    std::unordered_map<std::string, ExtModelArrayFlagged> store;

    template <const TextTag& fieldName, bool isReadWrite> friend class ModelArrayAccessor;
};

}

#endif /* MARSTORE_HPP */
