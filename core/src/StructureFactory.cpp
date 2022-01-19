/*!
 * @file StructureFactory.cpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/StructureFactory.hpp"

#include <ncFile.h>
#include <ncGroup.h>
#include <ncGroupAtt.h>
#include <stdexcept>
#include <string>

namespace Nextsim {

std::shared_ptr<IStructure> StructureFactory::generate(const std::string& structureName)
{
    ModuleLoader& loader = ModuleLoader::getLoader();
    std::string iStruct = "Nextsim::IStructure";
    std::shared_ptr<IStructure> shst;
    for (auto struc : loader.listImplementations(iStruct)) {
        loader.setImplementation(iStruct, struc);
        if (loader.getImplementation<IStructure>().structureType() == structureName) {
            shst = std::move(loader.getInstance<IStructure>());
            return shst;
        }
    }
    // If we reach here, then no valid handlers of the named structure were
    // found. Throw a invalid argument exception.
    std::string what = "StructureSelector: Invalid structure name (";
    what += structureName + ") provided.";
    throw std::invalid_argument(what);
}

std::shared_ptr<IStructure> StructureFactory::generateFromFile(const std::string& filePath)
{
    netCDF::NcFile ncf(filePath, netCDF::NcFile::read);
    netCDF::NcGroup metaGroup(ncf.getGroup(IStructure::metadataNodeName()));
    netCDF::NcGroupAtt att = metaGroup.getAtt(IStructure::typeNodeName());
    int len = att.getAttLength();
    // Initialize a std::string of len, filled with zeros
    std::string structureName(len, '\0');
    // &str[0] gives access to the buffer, guaranteed by C++11
    att.getValues(&structureName[0]);
    ncf.close();

    return generate(structureName);
}

} /* namespace Nextsim */
