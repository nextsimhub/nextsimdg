/*!
 * @file NextsimModule.hpp
 *
 * @date Sep 13, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef NEXTSIMMODULE_HPP
#define NEXTSIMMODULE_HPP

#include "include/ConfigurationHelp.hpp"
#include "include/ConfiguredModule.hpp"

// Define the types the module should use (HelpMap and Config)
namespace Module {
using OptionMap = std::list<Nextsim::ConfigurationHelp>;
using HelpMap = std::map<std::string, OptionMap>;

using Config = Nextsim::ConfiguredModule;
}

#include "include/Module.hpp"

#endif /* NEXTSIMMODULE_HPP */
