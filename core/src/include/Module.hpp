/*!
 * @file Module.hpp
 *
 * @date Feb 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODULE_HPP
#define MODULE_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/ConfiguredModule.hpp"

#include <stdexcept>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <string>

namespace Module {

template <typename I> std::unique_ptr<I> getInstance();

template <typename I> I& getImplementation();

template <typename I> void setImplementation(const std::string&);

template <typename Int, typename Imp> std::unique_ptr<Int> newImpl()
{
    return std::unique_ptr<Int>(new Imp);
}

template <typename I, typename M> I& getImplTemplate() { return M::getImplementation(); }

template <typename M> void setImplTemplate(const std::string& implName)
{
    M::setImplementation(implName);
}

template <typename I, typename M> std::unique_ptr<I> getInstTemplate() { return M::getInstance(); }

typedef std::list<Nextsim::ConfigurationHelp> OptionMap;
typedef std::map<std::string, OptionMap> HelpMap;
using ConfigType = Nextsim::ConfigurationHelp::ConfigType;

template <typename I> HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

template <typename I> class Module {
public:
    typedef std::function<std::unique_ptr<I>()> fn;
    typedef std::map<std::string, fn> map;

    static void setExternalImplementation(fn generator)
    {
        spf = generator;
        staticInstance = std::move(spf());
    }

    static void setImplementation(const std::string& implName)
    {
        // setExternalImplementation() holds the common functionality
        try {
            setExternalImplementation(functionMap.at(implName));
        } catch (const std::out_of_range& oor) {
            std::throw_with_nested(std::runtime_error(
                "No implementation named " + implName + " found for Module " + moduleName()));
        }
    }

    static std::unique_ptr<I> getInstance() { return spf(); }

    static I& getImplementation() { return *staticInstance; }

    static std::list<std::string> listImplementations()
    {
        std::list<std::string> keys;
        for (auto entry : functionMap) {
            keys.push_back(entry.first);
        }
        return keys;
    }

    static std::string implementation()
    {
        typedef std::unique_ptr<I>(fnType)();
        // The name is not cached, so find the function in the map which
        // corresponds to spf. The hairy code is derived from
        // https://stackoverflow.com/a/35920804
        for (auto entry : functionMap) {
            if (*entry.second.template target<fnType*>() == *spf.template target<fnType*>()) {
                return entry.first;
            }
        }
        return ""; // spf should always be an entry in functionMap
    }

    static std::string moduleName();

private:
    static fn spf;
    static std::unique_ptr<I> staticInstance;
    static map functionMap;
};

template <typename I, typename M> void addToConfiguredModules()
{
    Nextsim::ConfiguredModule::configureModule(
        Module<I>::moduleName(), M::setImplementation, M::implementation);
}

template <typename I> std::string implementation() { return Module<I>::implementation(); }

}
#endif /* MODULE_HPP */
