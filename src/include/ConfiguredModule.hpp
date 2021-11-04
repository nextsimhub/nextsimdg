/*
 * @file ConfiguredModule.hpp
 *
 * @date Oct 29, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_CONFIGUREDMODULE_HPP
#define SRC_INCLUDE_CONFIGUREDMODULE_HPP

namespace Nextsim {

class ConfiguredModule {
 public:
  ConfiguredModule();
  virtual ~ConfiguredModule();

  static void parseConfigurator();
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_CONFIGUREDMODULE_HPP */
