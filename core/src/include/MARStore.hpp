/*!
 * @file MARBackingStore.hpp
 *
 * @date 22 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MARSTORE_HPP
#define MARSTORE_HPP

#include "include/TextTag.hpp"

#include <list>
#include <string>
#include <unordered_map>

namespace Nextsim {

class ModelArray;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;

class MARStore {
public:
    void registerArray(const std::string& field, ModelArray* ptr, bool isReadWrite = false)
    {
        if (!isReadWrite) {
            storeRO[field] = ptr;
        } else {
            storeRW[field] = ptr;
            if (waitingRW.count(field)) {
                for (ModelArrayReference* ref : waitingRW.at(field)) {
                    *ref = ptr;
                }
            }
        }
        if (waitingRO.count(field)) {
            for (ModelArrayConstReference* ref : waitingRO.at(field)) {
                *ref = ptr;
            }
        }
    }
private:
    ModelArray* getFieldAddr(const std::string& field, ModelArrayReference& ptr)
    {
        // Add this address to the waiting list for RW fields.
        if (waitingRW.count(field))
            waitingRW.at(field).push_back(&ptr);
        else
            waitingRW[field] = { &ptr };
        return storeRW.count(field) ? ptr = storeRW.at(field) : nullptr;
    }

    const ModelArray* getFieldAddr(const std::string& field, ModelArrayConstReference& ptr)
    {
        // Add this address to the waiting list for RO fields.
        if (waitingRO.count(field))
            waitingRO.at(field).push_back(&ptr);
        else
            waitingRO[field] = { &ptr };

        if (storeRO.count(field)) {
            ptr = storeRO.at(field);
            return ptr;
        } else if (storeRW.count(field)) {
            ptr = storeRW.at(field);
            return ptr;
        } else {
            return nullptr;
        }
    }


    void removeReference(ModelArrayConstReference* ptr)
    {
        for (auto& [name, list] : waitingRO) {
            list.remove(ptr);
            if (list.size() == 0)
                waitingRO.erase(name);
        }
    }

    void removeReference(ModelArrayReference* ptr)
    {
        for (auto& [name, list] : waitingRW) {
            list.remove(ptr);
            if (list.size() == 0)
                waitingRW.erase(name);
        }
        for (auto& [name, list] : waitingRO) {
            list.remove(const_cast<ModelArrayConstReference*>(ptr));
            if (list.size() == 0)
                waitingRO.erase(name);
        }
    }

    std::unordered_map<std::string, ModelArray*> storeRO;
    std::unordered_map<std::string, ModelArray*> storeRW;
    std::unordered_map<std::string, std::list<ModelArrayReference*>> waitingRW;
    std::unordered_map<std::string, std::list<ModelArrayConstReference*>> waitingRO;

    template <const TextTag& fieldName, bool isReadWrite>
    friend class ModelArrayRef;
};

}

#endif /* MARSTORE_HPP */
