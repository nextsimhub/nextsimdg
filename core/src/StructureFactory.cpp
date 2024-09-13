/*!
 * @file StructureFactory.cpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#include "include/StructureFactory.hpp"

#include "include/IStructure.hpp"
#include "include/Module.hpp"

#include "include/RectGridIO.hpp"

#include "include/ParaGridIO.hpp"

#include <ncFile.h>
#include <ncGroup.h>
#include <ncGroupAtt.h>
#include <stdexcept>
#include <string>

namespace Nextsim {

std::string structureNameFromFile(const std::string& filePath)
{
    std::string structureName;

    try {
        netCDF::NcFile ncf(filePath, netCDF::NcFile::read);
        netCDF::NcGroup metaGroup(ncf.getGroup(IStructure::structureNodeName()));
        netCDF::NcGroupAtt att = metaGroup.getAtt(IStructure::typeNodeName());
        int len = att.getAttLength();
        // Initialize a std::string of len, filled with zeros
        structureName = std::string(len, '\0');
        // &str[0] gives access to the buffer, guaranteed by C++11
        att.getValues(&structureName[0]);
        ncf.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filePath;
        throw std::runtime_error(ncWhat);
    }

    return structureName;
}

#ifdef USE_MPI
ModelState StructureFactory::stateFromFile(const std::string& filePath, ModelMetadata& metadata)
#else
ModelState StructureFactory::stateFromFile(const std::string& filePath)
#endif
{
    std::string structureName = structureNameFromFile(filePath);
    // TODO There must be a better way
    if (RectangularGrid::structureName == structureName) {
        Module::setImplementation<IStructure>("Nextsim::RectangularGrid");
        RectangularGrid gridIn;
        gridIn.setIO(new RectGridIO(gridIn));
#ifdef USE_MPI
        return gridIn.getModelState(filePath, metadata);
#else
        return gridIn.getModelState(filePath);
#endif
    } else if (ParametricGrid::structureName == structureName) {
        Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
        ParametricGrid gridIn;
        gridIn.setIO(new ParaGridIO(gridIn));
#ifdef USE_MPI
        return gridIn.getModelState(filePath, metadata);
#else
        return gridIn.getModelState(filePath);
#endif
    } else {
        throw std::invalid_argument(
            std::string("fileFromName: structure not implemented: ") + structureName);
    }
    throw std::invalid_argument(std::string("fileFromName: structure not implemented: ")
        + structureName + "\nAlso, how did you get here?");
    return ModelState();
}

void StructureFactory::fileFromState(
    const ModelState& state, const ModelMetadata& meta, const std::string& filePath, bool isRestart)
{
    std::string structureName = Module::getImplementation<IStructure>().structureType();

    if (RectangularGrid::structureName == structureName) {
        RectangularGrid gridOut;
        gridOut.setIO(new RectGridIO(gridOut));
        gridOut.dumpModelState(state, meta, filePath, isRestart);
    } else if (ParametricGrid::structureName == structureName) {
        ParametricGrid gridOut;
        gridOut.setIO(new ParaGridIO(gridOut));
        gridOut.dumpModelState(state, meta, filePath, isRestart);
    } else {
        throw std::invalid_argument(
            std::string("fileFromName: structure not implemented: ") + structureName);
    }
}

} /* namespace Nextsim */
