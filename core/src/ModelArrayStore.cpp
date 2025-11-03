/*!
 * @author  Robert Jendersie <robert.jendersie@ovgu.de>
 */

#include "include/ModelArrayStore.hpp"

namespace Nextsim {

HostView ModelArrayStore::ExtModelArray::hostView()
{
    return makeKokkosHostView(modelArray.getDataRef());
}

ConstHostView ModelArrayStore::ExtModelArray::hostView() const
{
    return makeKokkosHostView(modelArray.data());
}

void ModelArrayStore::resize_arrays()
{
    was_resized = true;

    for (auto& [name, extArrFlagged] : store) {
        ExtModelArray& extArr = extArrFlagged.extModelArray;
        extArr.modelArray.resize();

        extArr.deviceView = makeKokkosDeviceView(name, extArr.modelArray.getDataRef());
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
    return store.try_emplace(field, RW, false).first->second.extModelArray;
}

const ModelArrayStore::ExtModelArray& ModelArrayStore::getRO(const std::string& field)
{
    return store.try_emplace(field, RO, false).first->second.extModelArray;
}
}