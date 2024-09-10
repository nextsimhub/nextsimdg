/*!
 * @file MARStore.hpp
 *
 * @date 30 Aug 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MARSTORE_HPP
#define MARSTORE_HPP

#include <list>
#include <string>
#include <unordered_map>

#include <iostream> // FIXME remove me

struct TextTag;

namespace Nextsim {

class ModelArray;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;

class ModelArrayReferenceStore {
public:
    ~ModelArrayReferenceStore() { std::cout << this << ":" << "~ModelArrayReferenceStore" << std::endl; }
    void registerArray(const std::string& field, ModelArray* ptr, bool isReadWrite = false)
    {
        std::cout << this << ":" << "MARS::registerArray(" << field << ")" << std::endl;
        // Clean the old array, if set
        if (storeRW.count(field)) {
            // If the pointer was in the RW array, null it
            storeRW[field] = nullptr;
            // Set any RW references with this field to nullptr
            auto range = referencesRW.equal_range(field);
            for (auto it = range.first; it != range.second; ++it) {
                *(it->second) = nullptr;
            }
        } else if (storeRO.count(field)) {
            storeRO[field] = nullptr;
        }
        // Set any RO references with this field to nullptr
        auto range = referencesRO.equal_range(field);
        for (auto it = range.first; it != range.second; ++it) {
            *(it->second) = nullptr;
        }

        // Set the new reference based on the passed access boolean
        if (!isReadWrite) {
            storeRO[field] = ptr;
        } else {
            storeRW[field] = ptr;
            auto range = referencesRW.equal_range(field);
            for (auto it = range.first; it != range.second; ++it) {
                *(it->second) = ptr;
            }
        }
        range = referencesRO.equal_range(field);
        for (auto it = range.first; it != range.second; ++it) {
            *(it->second) = ptr;
        }
    }

    std::unordered_map<std::string, const ModelArray*> getAllData() const
    {
        std::cout << this << ":" << "MARS::getAllData()" << std::endl;
        std::unordered_map<std::string, const ModelArray*> dataMap;

        for (auto entry : storeRW) {
            dataMap.insert(entry);
        }
        for (auto entry : storeRO) {
            dataMap.insert(entry);
        }
        return dataMap;
    }

private:
    ModelArray* getFieldAddr(const std::string& field, ModelArrayReference& ptr)
    {
        std::cout << this << ":" << "MARS::getFieldAddr(" << field << ")" << std::endl;
        // Add this address to the waiting list for RW fields.
        referencesRW.insert({ field, &ptr });
        return storeRW.count(field) ? ptr = storeRW.at(field) : nullptr;
    }

    const ModelArray* getFieldAddr(const std::string& field, ModelArrayConstReference& ptr)
    {
        std::cout << this << ":" << "MARS::getFieldAddr(" << field << ")" << std::endl;
        referencesRO.insert({ field, &ptr });
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

    void removeReference(const std::string& field, ModelArrayConstReference& ptr)
    {
        std::cout << this << ":" << "MARS::removeReference(" << field << ")" << std::endl;
        auto range = referencesRO.equal_range(field);
        for (auto it = range.first; it != range.second;)
            if (it->second == &ptr)
                it = referencesRO.erase(it);
            else
                ++it;
    }

    void removeReference(const std::string& field, ModelArrayReference& ptr)
    {
        std::cout << this << ":" << "MARS::getFieldAddr(" << field << ")" << std::endl;
        auto range = referencesRW.equal_range(field);
        for (auto it = range.first; it != range.second;)
            if (it->second == &ptr)
                it = referencesRW.erase(it);
            else
                ++it;
    }

    std::unordered_map<std::string, ModelArray*> storeRO;
    std::unordered_map<std::string, ModelArray*> storeRW;
    std::unordered_multimap<std::string, ModelArrayReference*> referencesRW;
    std::unordered_multimap<std::string, ModelArrayConstReference*> referencesRO;
    template <const TextTag& fieldName, bool isReadWrite> friend class ModelArrayRef;
};

}

#endif /* MARSTORE_HPP */
