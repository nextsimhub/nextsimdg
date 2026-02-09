/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/ModelArrayStore.hpp"
#include "include/Logged.hpp"
#include "include/ModelArrayAccessor.hpp"
#include <iostream>
namespace Nextsim {

#ifdef USE_KOKKOS

/*************************************************************/
HostViewMA ModelArrayStore::ExtModelArray::hostView()
{
    assert(modelArray.trueSize() > 0 && "ModelArray is allocated");
    return makeKokkosHostView(modelArray.getDataRef());
}

const DeviceViewMA& ModelArrayStore::ExtModelArray::deviceView()
{
    assert(modelArray.trueSize() > 0 && "ModelArray is allocated");
    if (!m_deviceModelArray.m_deviceView.is_allocated()) {
        m_deviceModelArray.m_deviceView = makeKokkosDeviceView(name, modelArray.getDataRef());
    }

    return m_deviceModelArray.m_deviceView;
}

DeviceModelArray& ModelArrayStore::ExtModelArray::deviceModelArray()
{
    // handles initialization
    deviceView();
    return m_deviceModelArray;
}
#endif

/*************************************************************/
ModelArrayStore::ExtModelArrayFlagged::ExtModelArrayFlagged(
    const std::string& name, bool _isReadWrite, bool _isRegistered)
    : isReadWrite(_isReadWrite)
    , isRegistered(_isRegistered)
{
    extModelArray.name = name;
}

/*************************************************************/
/*std::unordered_map<std::string, const ModelArray*> ModelArrayStore::getAllData() const
{
    std::unordered_map<std::string, const ModelArray*> dataMap;

    for (const auto& [name, extArrFlagged] : store) {
        dataMap.emplace(name, &extArrFlagged.extModelArray.modelArray);
    }

    return dataMap;
}

std::unordered_map<std::string, ModelArrayAccessorBase<RO>> ModelArrayStore::getAllData() const
{
    std::unordered_map<std::string, ModelArrayAccessorBase<RO>> dataMap;

    for (const auto& [name, extArrFlagged] : store) {
        // Internally ModelArrayAccessor always holds a mutable reference but the RO variant only
exposes
        // the data as const so this is safe.
        dataMap.emplace(name,
ModelArrayAccessorBase<RO>(const_cast<ExtModelArrayFlagged&>(extArrFlagged)));
    }

    return dataMap;
}*/

/*************************************************************/
bool ModelArrayStore::checkAllRegistered() const
{
    bool b = true;
    for (const auto& [_, extArrFlagged] : store) {
        if (!extArrFlagged.isRegistered) {
            Logged::warning("Field \"" + extArrFlagged.extModelArray.name
                + "\" was not properly initialized by a registering ModelArrayAccessor.\n");
            b = false;
        }
    }

    return b;
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