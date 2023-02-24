/*!
 * @file IMARBackingStore.hpp
 *
 * @date 22 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IMARBACKINGSTORE_HPP
#define IMARBACKINGSTORE_HPP

#include <string>
#include <unordered_map>
#include <list>

#include <iostream> // FIXME remove me

namespace Nextsim {

class ModelArray;
//class ModelArrayRef;

typedef ModelArray* ModelArrayReference;
typedef const ModelArray* ModelArrayConstReference;

class IMARBackingStore {
public:
    ModelArray* getFieldAddr(const std::string& field, ModelArrayReference& ptr)
    {
        std::cerr << "IMARBS::gFA non-const: ";
        try {
            std::cerr << "store = {";
            for (const auto& entry : storeRW)
                std::cerr << "{" << entry.first << "," << entry.second << "},";
            std::cerr << "}" << std::endl;
            ptr = storeRW.at(field);
            return ptr;
        } catch (std::out_of_range& oore) {
            // If there is no entry for the field key, then add it to the waiting list
            if (waitingRW.count(field))
                waitingRW.at(field).push_back(&ptr);
            else
                waitingRW[field] = {&ptr};
            std::cerr << " RW MARs waiting on " << field << " = " << waitingRW.at(field).size() << std::endl;
            return nullptr;
        }
    }

    const ModelArray* getFieldAddr(const std::string& field, ModelArrayConstReference& ptr)
    {
        std::cerr << "IMARBS::gFA const: ";
        try {
            std::cerr << "store = {";
            for (const auto& entry : storeRO)
                std::cerr << "{" << entry.first << "," << entry.second << "},";
            std::cerr << "}" << std::endl;
            ptr = storeRO.at(field);
            return ptr;
        } catch (std::out_of_range& oore_ro) {
            try {
                std::cerr << "storeRW = {";
                for (const auto& entry : storeRW)
                    std::cerr << "{" << entry.first << "," << entry.second << "},";
                std::cerr << "}" << std::endl;
                ptr = storeRW.at(field);
                return ptr;
            } catch (std::out_of_range& oore) {
                // If there is no entry for the field key, then add it to the waiting list
                if (waitingRO.count(field))
                    waitingRO.at(field).push_back(&ptr);
                else
                    waitingRO[field] = {&ptr};
                std::cerr << " RO MARs waiting on " << field << " = " << waitingRO.at(field).size() << std::endl;
                return nullptr;
            }
        }
    }

    void registerArray(const std::string& field, ModelArray* ptr, bool isReadWrite = false)
    {
        if (!isReadWrite) {
            storeRO[field] = ptr;
        } else {
            storeRW[field] = ptr;
            if (waitingRW.count(field)) {
                std::cerr << "RW waiting arrays on field " << field << " = " << waitingRW.at(field).size() << std::endl;
                for (ModelArrayReference* ref : waitingRW.at(field)) {
                    *ref = ptr;
                }
            }
        }
        if (waitingRO.count(field)) {
            std::cerr << "RO waiting arrays on field " << field << " = " << waitingRO.at(field).size() << std::endl;
            for (ModelArrayConstReference* ref : waitingRO.at(field)) {
                *ref = ptr;
            }
        }
    }

    void removeReference(ModelArrayConstReference* ptr)
    {
        std::cerr << "Removing const reference " << ptr << std::endl;
        for (auto& [name, list] : waitingRO) {
            list.remove(ptr);
            if (list.size() == 0)
                waitingRO.erase(name);
        }
    }

    void removeReference(ModelArrayReference *ptr)
    {
        std::cerr << "Removing non-const reference " << ptr << std::endl;
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

private:
    std::unordered_map<std::string, ModelArray*> storeRO;
    std::unordered_map<std::string, ModelArray*> storeRW;
    std::unordered_map<std::string, std::list<ModelArrayReference*>> waitingRW;
    std::unordered_map<std::string, std::list<ModelArrayConstReference*>> waitingRO;
};

}

#endif /* IMARBACKINGSTORE_HPP */
