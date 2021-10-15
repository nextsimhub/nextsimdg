/*!
 * @file Configured.hpp
 *
 * @date Oct 8, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONFIGURED_HPP
#define SRC_INCLUDE_CONFIGURED_HPP

#include <set>
#include <boost/program_options.hpp>

namespace Nextsim {
class Configurator;

//! A base class to provide configuration infrastructure to class that can be configured.
class Configured {
public:
    //! The constructor adds the instantiated object to the list of objects to
    // be configured.
    Configured();
    virtual ~Configured() = default;

    //! Execute the parse() function on all objects that inherit from this
    // class.
    static void configureAll();

protected:
    /*!
     * Parse the config data in Configurator.
     *
     * Each class should implement its own logic that applies the configuration
     * in an override of this function. Calling this function from this class
     * will perform the actual parse on the config sources.
     */
    virtual void parse();

    /*!
     * Add a object to those that should be configured.
     *
     * Any instances of classes that inherit from this class call this
     * automatically when the base class constructor executes.
     *
     * @param object a raw pointer to an instance of a class that is derived
     * from this, which will be configured when configureAll() is run.
     */
    static void addConfiguredObject(Configured* object);

    /*!
     * Add an option to be configured for this class.
     *
     * A simple interface that wraps the boost::program_options underlying the
     * configuration options.
     *
     * @param name the name of the configuration option.
     * @param defaultValue the value that should be used if this option is not
     *          found in the configuration sources.
     * @param helpText some descriptive text about the option.
     */
    template<typename T>
    void addOption(const std::string name, const T& defaultValue, const std::string& helpText)
    {
        opt.add_options()
                (name.c_str(), boost::program_options::value<T>()->default_value(defaultValue), helpText.c_str())
                ;
    }

    //! Retrieve the value of an option, as configured.
    template<typename T>
    T retrieveValue(const std::string& name)
    {
        return vm[name].as<T>();
    }

private:
    boost::program_options::options_description opt;
    boost::program_options::variables_map vm;
    // The derived objects might not be stored on the heap, so we cannot use
    // <memory> pointers to access them. Hence the use of raw pointers.
    static std::set<Configured*> configuredObjects;
};
}

#endif /* SRC_INCLUDE_CONFIGURED_HPP */
