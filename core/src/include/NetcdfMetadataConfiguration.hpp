/*!
 * @file NetcdfMetadataConfiguration.hpp
 *
 * @date Aug 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef NETCDFMETADATACONFIGURATION_HPP
#define NETCDFMETADATACONFIGURATION_HPP

#include "include/ConfigMap.hpp"
#include "include/Configurator.hpp"

namespace Nextsim {

/*!
 * An implementation of Configurator::AdditionalConfiguration which fetches
 * configuration from the /metadata/configuration node of the named netCDF file.
 */
class NetcdfMetadataConfiguration : public Configurator::AdditionalConfiguration {
public:
    std::stringstream read(const std::string& source) override;
};

} /* namespace Nextsim */

#endif /* NETCDFMETADATACONFIGURATION_HPP */
