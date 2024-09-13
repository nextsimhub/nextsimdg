/*!
 * @file Module.hpp
 *
 * @date Feb 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODULE_HPP
#define MODULE_HPP

#include <list>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

namespace Module {

/*!
 * A class that records the names of all Modules in this executable.
 */
class ImplementationNames {
public:
    using Container = std::map<std::string, std::string>;
    static const Container& getAll() { return getNames(); }
    static void set(const std::string& interface, const std::string& implementation)
    {
        getNames().insert({ interface, implementation });
    }
    static std::string get(const std::string& interface) { return getNames().at(interface); }

private:
    static Container& getNames()
    {
        static Container cache;
        return cache;
    }
};

template <typename I, typename C, typename H> class Module;

template <typename I, typename C, typename H> std::unique_ptr<I> getInstance() { return Module<I, C, H>::getInstance(); }

template <typename I, typename C, typename H> I& getImplementation() { return Module<I, C, H>::getImplementation(); }

template <typename I, typename C, typename H> void setImplementation(const std::string& impl) { Module<I, C, H>::setImplementation(impl); }

template <typename Int, typename Imp> std::unique_ptr<Int> newImpl()
{
    return std::unique_ptr<Int>(new Imp);
}

template <typename I, typename C, typename H> H& getHelpRecursive(H& map, bool getAll)
{
    return Module<I, C, H>::getHelpRecursive(map, getAll);
}

template <typename I, typename C, typename H> std::string implementation() { return Module<I, C, H>::implementation(); }

template <typename I, typename C, typename H> class Module {
public:
    using fn = std::unique_ptr<I> (*)();
    using map = std::map<std::string, fn>;

    static void setExternalImplementation(fn generator)
    {
        setExternalImplementation(generator, true);
    }

    static void setImplementation(const std::string& implName)
    {
        setImplementation(implName, true);
    }

    static std::unique_ptr<I> getInstance() { return getInstance(true); }

    static std::unique_ptr<I>& getUniqueInstance()
    {
        static std::unique_ptr<I> staticInstance = std::move(getInstance(false));
        return staticInstance;
    }

    static I& getImplementation()
    {
        if (!getUniqueInstance().get())
            throw std::logic_error("Bad reference implementation in " + moduleName());
        return *getUniqueInstance();
    }

    static std::list<std::string> listImplementations()
    {
        std::list<std::string> keys;
        for (const auto& entry : functionMap()) {
            keys.push_back(entry.first);
        }
        return keys;
    }

    static std::string implementation()
    {
        typedef std::unique_ptr<I>(fnType)();
        // The name is not cached, so find the function in the map which
        // corresponds to getGenerationFunction.
        for (const auto& entry : functionMap()) {
            if (entry.second == getGenerationFunction()) {
                return entry.first;
            }
        }
        throw std::out_of_range("Module<" + moduleName() + ">: implementation not found.");
        return ""; // getGenerationFunction should always be an entry in functionMap
    }

    static std::string moduleName();

    static H& getHelpRecursive(H& helpMap, bool getAll);

private:
    static fn& getGenerationFunction();
    static const map& functionMap();

    static bool& isConfigured()
    {
        static bool isConfiguredBool = false;
        return isConfiguredBool;
    }

    static std::unique_ptr<I> getInstance(bool setStaticInstance)
    {
        if (!isConfigured()) {
            isConfigured() = true;
            std::string implName = C::getImpl(moduleName());
            if (!implName.empty()) {
                setImplementation(implName, setStaticInstance);
            }
        }
        return getGenerationFunction()();
    }

    static void setImplementation(const std::string& implName, bool setStaticInstance)
    {
        // setExternalImplementation() holds the common functionality
        try {
            setExternalImplementation(functionMap().at(implName), setStaticInstance);
            ImplementationNames::set(moduleName(), implementation());
        } catch (const std::out_of_range& oor) {
            std::throw_with_nested(std::runtime_error(
                "No implementation named " + implName + " found for Module " + moduleName()));
        }
    }

    static void setExternalImplementation(fn generator, bool setStaticInstance)
    {
        getGenerationFunction() = generator;
        if (setStaticInstance)
            getUniqueInstance() = std::move(getGenerationFunction()());
    }
};

}
#endif /* MODULE_HPP */
