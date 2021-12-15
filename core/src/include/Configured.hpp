/*!
 * @file Configured.hpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONFIGURED_HPP
#define SRC_INCLUDE_CONFIGURED_HPP

#include "Configurator.hpp"

#include <boost/program_options.hpp>
#include <map>
#include <string>
#include <typeinfo>

namespace Nextsim {

//! A base class necessary for the tryConfigure() functions to work.
class ConfiguredBase {
public:
    virtual ~ConfiguredBase() = default;
    //! The configuration function.
    virtual void configure() = 0;
};

/*!
 * A base class to provide configuration infrastructure to classes that can be
 * configured.
 */
template <typename C> class Configured : public ConfiguredBase {
public:
    Configured() = default;
    virtual ~Configured() = default;

    //! The configuration function.
    virtual void configure() = 0;

    /*!
     * @brief Template function for conditionally configuring references.
     *
     * @detailed Pass any class to this function (or its pointer equivalent).
     * If it is a derived class of Configured, the overridden configure()
     * function will be called. If it is not, the function will do nothing,
     * gracefully.
     *
     * @param ref A reference to the class on which to attempt configuration.
     */
    template <typename T> static void tryConfigure(T& ref);

    /*!
     * @brief Template function for conditionally configuring classes via a
     * pointer.
     *
     * @detailed Pass any class to this function (or its reference equivalent).
     * If it is a derived class of Configured, the overridden configure()
     * function will be called. If it is not, the function will do nothing,
     * gracefully.
     *
     * @param ptr A pointer to the class on which to attempt configuration.
     */
    template <typename T> static void tryConfigure(T* ptr);

    /*!
     * @brief Gets the value of the configuration with a given name from the
     * default Configurator.
     *
     * @param name Name of the configuration option to fetch.
     * @param defaultValue Default value to apply if the configuration is not
     * found.
     */
    template <typename T>
    static inline T getConfiguration(const std::string& name, const T& defaultValue)
    {
        boost::program_options::options_description opt;
        addOption(name, defaultValue, opt);
        return retrieveValue<T>(name, opt);
    }

    //! Clear the configuration map. Usually used only in test suites.
    static void clearConfigurationMap() { singleOptions.clear(); }

    //! A per-class static map to provide compile-time checking of configuration keys.
    static const std::map<int, std::string> keyMap;

protected:
    /*!
     * @brief Adds an option to the per-class option map.
     *
     * @param name Name of the option to add.
     * @param defaultValue Default value to apply if the configuration is not
     * found.
     */
    template <typename T> void addOption(const std::string& name, const T& defaultValue)
    {
        addOption(name, defaultValue, singleOptions[name]);
    }

    /*!
     * @brief Retrieves a configured value of a single option.
     *
     * @param name Name of the configuration option to fetch.
     */
    template <typename T> T retrieveValue(const std::string& name)
    {
        return retrieveValue<T>(name, singleOptions.at(name));
    }

private:
    template <typename T>
    static void addOption(const std::string& name, const T& defaultValue,
        boost::program_options::options_description& opt)
    {
        opt.add_options()(
            name.c_str(), boost::program_options::value<T>()->default_value(defaultValue), "");
    }

    template <typename T>
    static T retrieveValue(
        const std::string& name, boost::program_options::options_description& opt)
    {
        return Configurator::parse(opt)[name].as<T>();
    }

    static std::map<std::string, boost::program_options::options_description> singleOptions;
};

template <typename C>
std::map<std::string, boost::program_options::options_description> Configured<C>::singleOptions;

/*!
 * @brief Template function for conditionally configuring references.
 *
 * @detailed Pass any class to this function (or its pointer equivalent). If it
 * is a derived class of Configured, the overridden configure() function will
 * be called. If it is not, the function will do nothing, gracefully.
 *
 * @param ref A reference to the class on which to attempt configuration.
 */
template <typename C> template <typename T> void Configured<C>::tryConfigure(T& ref)
{
    try {
        dynamic_cast<ConfiguredBase&>(ref).configure();
    } catch (const std::bad_cast& bc) {
        // Do nothing. If the reference is not a derived class of Configured, ignore it.
    }
}

/*!
 * @brief Template function for conditionally configuring classes via a
 * pointer.
 *
 * @detailed Pass any class to this function (or its reference equivalent). If
 * it is a derived class of Configured, the overridden configure() function
 * will be called. If it is not, the function will do nothing, gracefully.
 *
 * @param ptr A pointer to the class on which to attempt configuration.
 */
template <typename C> template <typename T> void Configured<C>::tryConfigure(T* ptr)
{
    ConfiguredBase* cfg = dynamic_cast<ConfiguredBase*>(ptr);
    if (cfg)
        cfg->configure();
}

/*!
 * @brief Template function for conditionally configuring references.
 *
 * @detailed Pass any class to this function (or its pointer equivalent). If it
 * is a derived class of Configured, the overridden configure() function will
 * be called. If it is not, the function will do nothing, gracefully.
 *
 * @param ref A reference to the class on which to attempt configuration.
 */
template <typename T> void tryConfigure(T& t) { Configured<int>::tryConfigure(t); }

/*!
 * @brief Template function for conditionally configuring classes via a
 * pointer.
 *
 * @detailed Pass any class to this function (or its reference equivalent).
 * If it is a derived class of Configured, the overridden configure()
 * function will be called. If it is not, the function will do nothing,
 * gracefully.
 *
 * @param ptr A pointer to the class on which to attempt configuration.
 */
template <typename T> void tryConfigure(T* p_t) { Configured<int>::tryConfigure(p_t); }
}

#endif /* SRC_INCLUDE_CONFIGURED_HPP */
