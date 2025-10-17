/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef COMMONRESTARTMETADATA_HPP
#define COMMONRESTARTMETADATA_HPP

#include "include/IStructure.hpp"

#include <ncFile.h>

namespace Nextsim {

class CommonRestartMetadata {
public:
    //! Writes the structure type to the root of the restart file for future
    //! retrieval.
    static netCDF::NcFile& writeStructureType(netCDF::NcFile& ncFile);
    //! Writes the standard restart file metadata to a metadata node.
    static netCDF::NcFile& writeRestartMetadata(netCDF::NcFile& ncFile);

    static const std::string timeNodeName() { return "time"; }

    static const std::string formattedName() { return "formatted"; }

    static const std::string unformattedName() { return "time_meta"; }

    static const std::string configurationNode() { return "configuration"; }

    static const std::string startTimeName() { return "start"; }

    static const std::string stopTimeName() { return "stop"; }

    static const std::string stepLengthName() { return "step"; }

    static const std::string runLengthName() { return "run_length"; }

    static const std::string missingDataName() { return "missing_value"; }
};

} /* namespace Nextsim */

#endif /* COMMONRESTARTMETADATA_HPP */
