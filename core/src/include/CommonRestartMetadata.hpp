/*!
 * @file CommonRestartMetadata.hpp
 *
 * @date Jun 30, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef COMMONRESTARTMETADATA_HPP
#define COMMONRESTARTMETADATA_HPP

#include "include/IStructure.hpp"
#include "include/ModelMetadata.hpp"

#include <ncFile.h>
#include <ncGroup.h>

namespace Nextsim {

class CommonRestartMetadata {
public:
    //! Writes the structure type to the root of the restart file for future
    //! retrieval.
    static netCDF::NcGroup& writeStructureType(
        netCDF::NcFile& rootGroup, const ModelMetadata& metadata);
    //! Writes the standard restart file metadata to a metadata node.
    static netCDF::NcGroup& writeRestartMetadata(
        netCDF::NcGroup& metaGroup, const ModelMetadata& metadata);

    static const std::string timeNodeName() { return "time"; }

    static const std::string formattedName() { return "formatted"; }

    static const std::string unformattedName() { return "time"; }

    static const std::string configurationNode() { return "configuration"; }
};

} /* namespace Nextsim */

#endif /* COMMONRESTARTMETADATA_HPP */
