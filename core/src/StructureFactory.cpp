/*!
 * @file StructureFactory.cpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/StructureFactory.hpp"

#include "include/IStructureModule.hpp"

#include "include/DevGrid.hpp"
#include "include/DevGridIO.hpp"

#include "include/RectGridIO.hpp"
#include "include/RectangularGrid.hpp"

#include <ncFile.h>
#include <ncGroup.h>
#include <ncGroupAtt.h>
#include <stdexcept>
#include <string>

namespace Nextsim {

//std::shared_ptr<IStructure> StructureFactory::generate(const std::string& structureName)
//{
//    std::string iStruct = "Nextsim::IStructure";
//    std::shared_ptr<IStructure> shst;
//    for (auto struc : Module::IStructureModule::listImplementations()) {
//        Module::setImplementation<IStructure>(struc);
//        if (Module::getImplementation<IStructure>().structureType() == structureName) {
//            shst = std::move(Module::getInstance<IStructure>());
//
//            // TODO There must be a better way
//            if (shst->structureTypeCheck(DevGrid::structureName)) {
//                std::shared_ptr<DevGrid> shdg = std::dynamic_pointer_cast<DevGrid>(shst);
//                shdg->setIO(new DevGridIO(*shdg));
//            } else if (shst->structureTypeCheck(RectangularGrid::structureName)) {
//                std::shared_ptr<RectangularGrid> shrg
//                    = std::dynamic_pointer_cast<RectangularGrid>(shst);
//                shrg->setIO(new RectGridIO(*shrg));
//            }
//
//            return shst;
//        }
//    }
//    // If we reach here, then no valid handlers of the named structure were
//    // found. Throw a invalid argument exception.
//    std::string what = "StructureSelector: Invalid structure name (";
//    what += structureName + ") provided.";
//    throw std::invalid_argument(what);
//}

std::string structureNameFromFile(const std::string& filePath)
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

    return structureName;
}

//std::shared_ptr<IStructure> StructureFactory::generateFromFile(const std::string& filePath)
//{
//    return generate(structureNameFromFile(filePath));
//}

ModelState StructureFactory::stateFromFile(const std::string& filePath)
{
    std::string structureName = structureNameFromFile(filePath);
    for (auto struc : Module::IStructureModule::listImplementations()) {
        Module::setImplementation<IStructure>(struc);
        if (Module::getImplementation<IStructure>().structureType() == structureName) {
            // TODO There must be a better way
            if (DevGrid::structureName == structureName) {
                DevGrid gridIn;
                gridIn.setIO(new DevGridIO(gridIn));
                return gridIn.getModelState(filePath);
            } else if (RectangularGrid::structureName == structureName) {
                RectangularGrid gridIn;
                gridIn.setIO(new RectGridIO(gridIn));
                // return gridIn.getModelState(filePath);
                return ModelState();
            }
        }
    }
    // TODO: throw some kind of exception if we get here.
    return ModelState();
}

} /* namespace Nextsim */
