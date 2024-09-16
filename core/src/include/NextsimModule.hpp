/*!
 * @file NextsimModule.hpp
 *
 * @date Sep 13, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef NEXTSIMMODULE_HPP
#define NEXTSIMMODULE_HPP

#include "include/Module.hpp"
#include "include/ConfigurationHelp.hpp"
#include "include/ConfiguredModule.hpp"

namespace Nextsim {

using OptionMap = std::list<Nextsim::ConfigurationHelp>;
using HelpMap = std::map<std::string, OptionMap>;

template <typename I> class Module : public Modules::Module<I, ConfiguredModule, HelpMap> { };

template <typename I> std::unique_ptr<I> getInstance();
template <typename I> I& getImplementation();
template <typename I> void setImplementation(const std::string& impl);
template <typename I> HelpMap& getHelpRecursive(HelpMap& map, bool getAll);
template <typename I> std::string implementation();

}

#endif /* NEXTSIMMODULE_HPP */
