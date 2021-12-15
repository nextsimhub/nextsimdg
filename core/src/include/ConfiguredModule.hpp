/*!
 * @file ConfiguredModule.hpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONFIGUREDMODULE_HPP
#define SRC_INCLUDE_CONFIGUREDMODULE_HPP

#include <string>

namespace Nextsim {

class ConfiguredModule {
public:
    ConfiguredModule();
    virtual ~ConfiguredModule();

    static void parseConfigurator();

    static std::string addPrefix(const std::string& moduleName);

    static const std::string MODULE_PREFIX;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_CONFIGUREDMODULE_HPP */
