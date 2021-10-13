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
    Configured();
    virtual ~Configured() = default;

    static void configureAll();

protected:
    virtual void parse();

    static void addConfiguredObject(Configured*);

//    template<typename T>
//    void addOption(std::string& name, T& defaultValue, std::string& helpText)
//    {
//        opt.add_options()
//                (name.c_str(), boost::program_options::value<T>()->default_value(defaultValue), helpText.c_str())
//                ;
//    }

    template<typename T>
    void addOption(const std::string name, const T& defaultValue, const std::string& helpText)
    {
        opt.add_options()
                (name.c_str(), boost::program_options::value<T>()->default_value(defaultValue), helpText.c_str())
                ;
    }

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
