/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/ModelArrayStore.hpp"

namespace Nextsim {

#ifdef USE_KOKKOS

HostView ModelArrayStore::ExtModelArray::hostView()
{
    return makeKokkosHostView(modelArray.getDataRef());
}

DeviceView ModelArrayStore::ExtModelArray::deviceView()
{
    if (!m_deviceView.is_allocated()) {
        m_deviceView = makeKokkosDeviceView(name, modelArray.getDataRef());
    }

    return m_deviceView;
}

ModelArrayStore::ExtModelArrayFlagged::ExtModelArrayFlagged(
    const std::string& name, bool _isReadWrite, bool _isRegistered)
    : isReadWrite(_isReadWrite)
    , isRegistered(_isRegistered)
{
    extModelArray.name = name;
}

#endif

void ModelArrayStore::resize_arrays()
{
    was_resized = true;

    for (auto& [name, extArrFlagged] : store) {
        ExtModelArray& extArr = extArrFlagged.extModelArray;
        extArr.modelArray.resize();
    }
}

ModelArrayStore::ExtModelArray& ModelArrayStore::getRW(const std::string& field)
{
    if (auto it = store.find(field); it != store.end()) {
        ExtModelArrayFlagged& extArrFlagged = it->second;
        if (!extArrFlagged.isReadWrite) {
            if (extArrFlagged.isRegistered) {
                throw std::logic_error(
                    "Trying to access the read-only ModelArray \"" + field + "\" as read-write.");
            } else {
                // promote to read-write because we don't now the true restriction yet
                extArrFlagged.isReadWrite = true;
            }
        }
        return extArrFlagged.extModelArray;
    }

    // Regular emplace would be fine here since we know that it does not exist but try_emplace
    // has a more ergonomic signature.
    return store.try_emplace(field, field, RW, false).first->second.extModelArray;
}

ModelArrayStore::ExtModelArray& ModelArrayStore::getRO(const std::string& field)
{
    return store.try_emplace(field, field, RO, false).first->second.extModelArray;
}
}