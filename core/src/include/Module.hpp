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
    /*!
     * Returns a map between module name and current implementation name for all modules.
     */
    static Container getAll()
    {
        Container nameMap;
        for (auto entry : getNames()) {
            nameMap[entry.first()] = entry.second();
        }
        return nameMap;
    }

    /*!
     * Returns the implementation name of the named interface, else throws a std::out_of_range
     * exception.
     */
    static std::string get(const std::string& interface)
    {
        for (auto entry : getNames()) {
            if (entry.first() == interface) {
                return entry.second();
            }
        }
        throw std::out_of_range("No module found with name " + interface);
    }

private:
    using fn = std::string (*)();
    using InternalContainer = std::list<std::pair<fn, fn>>;
    static InternalContainer& getNames()
    {
        static InternalContainer cache;
        return cache;
    }

    static void set(fn intFunctionPtr, fn impFunctionPtr)
    {
        getNames().push_back({ intFunctionPtr, impFunctionPtr });
    }

    // Make all Module<>s a friend
    template <class Int> friend class Module;
};

template <typename Int, typename Imp> std::unique_ptr<Int> newImpl()
{
    return std::unique_ptr<Int>(new Imp);
}

template <typename I> class Module {
public:
    using fn = std::unique_ptr<I> (*)();
    using map = std::map<std::string, fn>;

    static void setExternalImplementation(fn generator)
    {
        setExternalImplementationInternal(generator, true);
    }

    static void setImplementation(const std::string& implName)
    {
        setImplementationInternal(implName, true);
    }

    static std::unique_ptr<I> getInstance() { return getInstanceInternal(true); }

    static std::unique_ptr<I>& getUniqueInstance()
    {
        static std::unique_ptr<I> staticInstance = std::move(getInstanceInternal(false));
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

    static HelpMap& getHelpRecursive(HelpMap& helpMap, bool getAll);

private:
    static std::string getDefaultImplementationName();
    static fn& getGenerationFunction()
    {
        static fn ptr = functionMap().at(getDefaultImplementationName());
        return ptr;
    }

    static const map& functionMap();

    static bool& isConfigured()
    {
        static bool isConfiguredBool = false;
        return isConfiguredBool;
    }

    static std::unique_ptr<I> getInstanceInternal(bool setStaticInstance)
    {
        if (!isConfigured()) {
            isConfigured() = true;
            std::string implName = Config::getImpl(moduleName());
            if (!implName.empty()) {
                setImplementationInternal(implName, setStaticInstance);
            }
        }
        return getGenerationFunction()();
    }

    static void setImplementationInternal(const std::string& implName, bool setStaticInstance)
    {
        // setExternalImplementation() holds the common functionality
        try {
            setExternalImplementationInternal(functionMap().at(implName), setStaticInstance);
            ImplementationNames::set(moduleName, implementation);
        } catch (const std::out_of_range& oor) {
            std::throw_with_nested(std::runtime_error(
                "No implementation named " + implName + " found for Module " + moduleName()));
        }
    }

    static void setExternalImplementationInternal(fn generator, bool setStaticInstance)
    {
        getGenerationFunction() = generator;
        if (setStaticInstance)
            getUniqueInstance() = std::move(getGenerationFunction()());
    }
};

template <typename I> std::unique_ptr<I> getInstance();
template <typename I> I& getImplementation();
template <typename I> void setImplementation(const std::string& impl);
template <typename I> HelpMap& getHelpRecursive(HelpMap& map, bool getAll);
template <typename I> std::string implementation();

}
#endif /* MODULE_HPP */
