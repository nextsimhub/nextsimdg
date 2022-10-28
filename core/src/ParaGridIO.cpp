/*!
 * @file ParaGridIO.cpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ParaGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/MissingData.hpp"
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

#include <map>
#include <string>

namespace Nextsim {

const std::map<std::string, ModelArray::Type> ParaGridIO::dimensionKeys = {
    { "xy", ModelArray::Type::H },
    { "xyz", ModelArray::Type::Z },

};

// Which dimensions are DG dimension, which could be legitimately missing
const std::map<ModelArray::Dimension, bool> ParaGridIO::isDG = {
    { ModelArray::Dimension::X, false },
    { ModelArray::Dimension::Y, false },
    { ModelArray::Dimension::Z, false },
    { ModelArray::Dimension::XCG, true },
    { ModelArray::Dimension::YCG, true },
    { ModelArray::Dimension::DG, true },
    { ModelArray::Dimension::DGSTRESS, true },
};

const std::map<ModelArray::Dimension, ModelArray::Type> ParaGridIO::dimCompMap = {
    { ModelArray::Dimension::DG, ModelArray::Type::DG },
    { ModelArray::Dimension::DGSTRESS, ModelArray::Type::DGSTRESS },
};
ModelState ParaGridIO::getModelState(const std::string& filePath)
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

    // Dimensions and DG components
    std::multimap<std::string, netCDF::NcDim> dimMap = dataGroup.getDims();

    for (auto entry : ModelArray::definedDimensions) {
        ModelArray::DimensionSpec& dimensionSpec = entry.second;
        netCDF::NcDim dim = dataGroup.getDim(dimensionSpec.name);
        if (dim.isNull()) {
            // Only throw if a finite volume dimension (X, Y, Z) is missing
            if (!isDG.at(entry.first)) {
                throw std::out_of_range(std::string("No dimension named ") + entry.second.name
                    + " found in " + filePath);
            }
        } else {
            if (dimCompMap.count(entry.first)) {
                ModelArray::setNComponents(dimCompMap.at(entry.first), dim.getSize());
            } else {
                dimensionSpec.length = dim.getSize();
            }
        }
    }

    // Get all vars in the data group, and load them into a new ModelState
    ModelState state;

    for (auto entry : dataGroup.getVars()) {
        const std::string& varName = entry.first;
        netCDF::NcVar& var = entry.second;
        // Determine the type from the dimensions
        std::vector<netCDF::NcDim> varDims = var.getDims();
        std::string dimKey = "";
        for (netCDF::NcDim& dim : varDims) {
            dimKey += dim.getName();
        }
        if (!dimensionKeys.count(dimKey)) {
            throw std::out_of_range(
                std::string("No ModelArray::Type corresponds to the dimensional key ") + dimKey);
        }
        state.data[varName](dimensionKeys.at(dimKey));
        ModelArray& data = state.data.at(varName);

        var.getVar(&data[0]);
    }

    ncFile.close();
    return state;
}

void ParaGridIO::dumpModelState(const ModelState& state, const ModelMetadata& metadata,
    const std::string& filePath, bool isRestart) const
{

    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);

    CommonRestartMetadata::writeStructureType(ncFile, metadata);
    netCDF::NcGroup metaGroup = ncFile.addGroup(IStructure::metadataNodeName());
    netCDF::NcGroup dataGroup = ncFile.addGroup(IStructure::dataNodeName());

    CommonRestartMetadata::writeRestartMetadata(metaGroup, metadata);

    // TODO dump grid data
    // Dump the dimensions and number of components
    std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
    for (auto entry : ModelArray::definedDimensions) {
        ModelArray::Dimension dim = entry.first;
        size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim)) : dimSz = entry.second.length;
        ncFromMAMap[dim] = dataGroup.addDim(entry.second.name, dimSz);
        // TODO Do I need to add data, even if it is just integers 0...n-1?
    }
    // Also create the sets of dimensions to be connected to the data fields
    std::map<ModelArray::Type, std::vector<netCDF::NcDim>> dimMap;
    for (auto entry : ModelArray::typeDimensions) {
        ModelArray::Type type = entry.first;
        std::vector<netCDF::NcDim> ncDims;
        for (ModelArray::Dimension& maDim : entry.second) {
            ncDims.push_back(ncFromMAMap.at(maDim));
        }
        dimMap[type] = ncDims;
    }

    std::set<std::string> restartFields
        = { hiceName, ciceName, hsnowName, ticeName }; // TODO and others
    // Loop through either the above list (isRestart) or all provided fields(!isRestart)
    for (auto entry : state.data) {
        if (!isRestart || restartFields.count(entry.first)) {
            // Get the type, then relevant vector of NetCDF dimensions
            ModelArray::Type type = entry.second.getType();
            std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
            netCDF::NcVar var(dataGroup.addVar(entry.first, netCDF::ncDouble, ncDims));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
            var.putVar(entry.second.getData());
        }
    }

    ncFile.close();
}
} /* namespace Nextsim */
