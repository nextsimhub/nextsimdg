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

const DeviceView& ModelArrayStore::ExtModelArray::deviceView()
{
    assert(modelArray.trueSize() > 0 && "ModelArray is allocated");
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