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

//! A base class to provide configuration infrastructure to class that can be configured.
class Configured {
public:
    Configured() = default;
    virtual ~Configured() = default;

    //! The configuration function.
    virtual void configure() = 0;

    //! Template function for conditionally configuring references.
    template <typename T> static void tryConfigure(T& ref)
    {
        try {
            dynamic_cast<Configured&>(ref).configure();
        } catch (const std::bad_cast& bc) {
            // Do nothing. If the reference is not a derived class of Configured, ignore it.
        }
    }

    //! Template function for conditionally configuring pointers.
    template <typename T> static void tryConfigure(T* ptr)
    {
        Configured* cfg = dynamic_cast<Configured*>(ptr);
        if (cfg)
            cfg->configure();
    }

    template <typename T>
    static inline T getConfiguration(const std::string& name, const T& defaultValue)
    {
        boost::program_options::options_description& opt;
        addOption(name, defaultValue, opt);
        return retrieveValue<T>(name, opt);
    }

protected:
    template <typename T> void addOption(const std::string& name, const T& defaultValue)
    {
        addOption(name, defaultValue, singleOptions[name]);
    }

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

    std::map<std::string, boost::program_options::options_description> singleOptions;
};
}

#endif /* SRC_INCLUDE_CONFIGURED_HPP */
