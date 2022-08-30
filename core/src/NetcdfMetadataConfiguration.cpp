/*!
 * @file NetcdfMetadataConfiguration.cpp
 *
 * @date Aug 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/NetcdfMetadataConfiguration.hpp"

#include <map>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

namespace Nextsim {

std::stringstream NetcdfMetadataConfiguration::read(const std::string& source)
{
    // Parse the current configuration to look for fileOptionName
    boost::program_options::options_description opt;
    opt.add_options()(
        source.c_str(), boost::program_options::value<std::string>()->default_value(""), "");
    std::string fileName = Configurator::parse(opt)[source].as<std::string>();

    netCDF::NcFile ncFile(fileName, netCDF::NcFile::read);
    netCDF::NcGroup configurationGroup(ncFile.getGroup("metadata").getGroup("configuration"));

    std::stringstream config;

    std::multimap<std::string, netCDF::NcVar> configs = configurationGroup.getVars();
    for (auto entry : configs) {
        config << entry.first << " = ";
        if (entry.second.getType() == netCDF::ncDouble) {
            double value;
            entry.second.getVar(&value);
            config << value;
        } else if (entry.second.getType() == netCDF::ncInt) {
            int value;
            entry.second.getVar(&value);
            config << value;
        } else if (entry.second.getType() == netCDF::ncUint) {
            unsigned value;
            entry.second.getVar(&value);
            config << value;
        } else if (entry.second.getType() == netCDF::ncString) {
            char* cValue;
            entry.second.getVar(&cValue);
            config << cValue;
        }
        config << std::endl;
    }

    ncFile.close();

    return config;
}

} /* namespace Nextsim */
