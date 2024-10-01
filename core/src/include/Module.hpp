/*!
 * @file Module.hpp
 *
 * @date 23 Sep 2024
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

    /*!
     * Sets the function to generate new instances of an implementation not included in the
     * relevant module.cfg file.
     *
     * @param generator A pointer to a function that returns a std::unique_ptr to an instance of
     * the template class.
     */
    static void setExternalImplementation(fn generator)
    {
        setExternalImplementationInternal(generator, true);
    }

    /*!
     * Sets the implementation of this module to one of the named implementations in the module.cfg
     * file.
     *
     * @param implName A string containing a name matching one of the implementations defined in
     * module.cfg.
     */
    static void setImplementation(const std::string& implName)
    {
        setImplementationInternal(implName, true);
    }

    /*!
     * Returns a std::unique_ptr to a new instance of the current implementation.
     */
    static std::unique_ptr<I> getInstance() { return getInstanceInternal(true); }

    /*!
     * Returns a std::unique_ptr to the static instance of the current implementation.
     */
    static std::unique_ptr<I>& getUniqueInstance(bool suppressInit = false)
    {
        static std::unique_ptr<I> staticInstance
            = std::move((suppressInit) ? std::unique_ptr<I>(nullptr) : getInstanceInternal(false));
        return staticInstance;
    }

    /*!
     * Returns a reference to the static instance of the current implementation.
     */
    static I& getImplementation()
    {
        if (!getUniqueInstance().get())
            throw std::logic_error("Bad reference implementation in " + moduleName());
        return *getUniqueInstance();
    }

    /*!
     * Returns a list of all the implementations named in the relevant module.cfg file.
     */
    static std::list<std::string> listImplementations()
    {
        std::list<std::string> keys;
        for (const auto& entry : functionMap()) {
            keys.push_back(entry.first);
        }
        return keys;
    }

    /*!
     * Returns the name of the current implementation , if it is listed in the relevant module.cfg
     * file. Otherwise returns an empty string.
     */
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
        /*
         *  If the generation function is not found in the function map, assume an external
         *  implementation and return an empty string.
         */
        return "";
    }

    //! Returns a string containing the name of the module.
    static std::string moduleName();

    //! Adds help information to the argument help map, according to the boolean getAll argument.
    //! Implementation dependent.
    static HelpMap& getHelpRecursive(HelpMap& helpMap, bool getAll);

    /*!
     * Finalizes the Module by setting both pointers to nullptr.
     */
    static void finalize()
    {
        getUniqueInstance(true) = nullptr;
        getGenerationFunction() = nullptr;
    }

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
            getUniqueInstance(true) = std::move(getGenerationFunction()());
    }
};

template <typename I> std::unique_ptr<I> getInstance();
template <typename I> I& getImplementation();
template <typename I> void setImplementation(const std::string& impl);
template <typename I> HelpMap& getHelpRecursive(HelpMap& map, bool getAll);
template <typename I> std::string implementation();
template <typename I> void finalize();

}
#endif /* MODULE_HPP */
