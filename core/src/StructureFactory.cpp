/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Kacper Kornet <kk562@cam.ac.uk>
 */

#include "include/StructureFactory.hpp"

#include "include/Finalizer.hpp"
#include "include/IStructure.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/RectGridIO.hpp"
#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif

#include <ncFile.h>
#include <ncGroupAtt.h>
#include <stdexcept>
#include <string>

namespace Nextsim {

std::string structureNameFromFile(const std::string& filePath)
{
    std::string structureName;

    try {
        netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
        netCDF::NcGroupAtt att = ncFile.getAtt(IStructure::structureNodeName());
        int len = att.getAttLength();
        // Initialize a std::string of len, filled with zeros
        structureName = std::string(len, '\0');
        // &str[0] gives access to the buffer, guaranteed by C++11
        att.getValues(&structureName[0]);
        ncFile.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filePath;
        throw std::runtime_error(ncWhat);
    }

    return structureName;
}

ModelState StructureFactory::stateFromFile(const std::string& filePath)
{
    Finalizer::registerUnique(Module::finalize<IStructure>);

#ifndef USE_XIOS
    std::string structureName = structureNameFromFile(filePath);
    // TODO There must be a better way
    if (RectangularGrid::structureName == structureName) {
        Module::setImplementation<IStructure>("Nextsim::RectangularGrid");
        RectangularGrid gridIn;
        gridIn.setIO(new RectGridIO(gridIn));
        return gridIn.getModelState(filePath);
    } else if (ParametricGrid::structureName == structureName) {
#endif
        Module::setImplementation<IStructure>("Nextsim::ParametricGrid");
        ParametricGrid gridIn;
        gridIn.setIO(new ParaGridIO(gridIn));
#ifdef USE_XIOS
        Xios& xiosHandler = Xios::getInstance();
        xiosHandler.close_context_definition();
#endif
        return gridIn.getModelState(filePath);
#ifndef USE_XIOS
    } else {
        throw std::invalid_argument(
            std::string("fileFromName: structure not implemented: ") + structureName);
    }
    throw std::invalid_argument(std::string("fileFromName: structure not implemented: ")
        + structureName + "\nAlso, how did you get here?");
#endif
    return ModelState();
}

void StructureFactory::fileFromState(
    const ModelState& state, const std::string& filePath, bool isRestart)
{
    std::string structureName = Module::getImplementation<IStructure>().structureType();

    if (RectangularGrid::structureName == structureName) {
        RectangularGrid gridOut;
        gridOut.setIO(new RectGridIO(gridOut));
        gridOut.dumpModelState(state, filePath, isRestart);
    } else if (ParametricGrid::structureName == structureName) {
        ParametricGrid gridOut;
        gridOut.setIO(new ParaGridIO(gridOut));
        gridOut.dumpModelState(state, filePath, isRestart);
    } else {
        throw std::invalid_argument(
            std::string("fileFromName: structure not implemented: ") + structureName);
    }
}

} /* namespace Nextsim */
