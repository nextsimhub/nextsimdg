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
    void parse();

    static void addConfiguredObject(Configured*);

    boost::program_options::options_description opt;
    boost::program_options::variables_map vm;

private:
    // The derived objects might not be stored on the heap, so we cannot use
    // <memory> pointers to access them. Hence the use of raw pointers.
    static std::set<Configured*> configuredObjects;
};
}

#endif /* SRC_INCLUDE_CONFIGURED_HPP */
