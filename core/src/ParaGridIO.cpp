/*!
 * @file ParaGridIO.cpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ParaGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/MissingData.hpp"
#include "include/NZLevels.hpp"
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

#include <algorithm>
#include <map>
#include <string>

namespace Nextsim {

const std::map<std::string, ModelArray::Type> ParaGridIO::dimensionKeys = {
    { "yx", ModelArray::Type::H },
    { "zyx", ModelArray::Type::Z },
    { "yxdg_comp", ModelArray::Type::DG },
    { "yxdgstress_comp", ModelArray::Type::DGSTRESS },
    { "ycgxcg", ModelArray::Type::CG },
    { "yvertexxvertexncoords", ModelArray::Type::VERTEX },
};

// Which dimensions are DG dimension, which could be legitimately missing
const std::map<ModelArray::Dimension, bool> ParaGridIO::isDG = {
    // clang-format off
    { ModelArray::Dimension::X, false },
    { ModelArray::Dimension::Y, false },
    { ModelArray::Dimension::Z, false },
    { ModelArray::Dimension::XCG, true },
    { ModelArray::Dimension::YCG, true },
    { ModelArray::Dimension::DG, true },
    { ModelArray::Dimension::DGSTRESS, true },
    // NCOORDS is a number of components, but not in the same way as the DG components.
    { ModelArray::Dimension::NCOORDS, false },
    // clang-format on
};

std::map<ModelArray::Dimension, ModelArray::Type> ParaGridIO::dimCompMap;

void ParaGridIO::makeDimCompMap()
{
    dimCompMap = {
        { ModelArray::componentMap.at(ModelArray::Type::DG), ModelArray::Type::DG },
        { ModelArray::componentMap.at(ModelArray::Type::DGSTRESS), ModelArray::Type::DGSTRESS },
        { ModelArray::componentMap.at(ModelArray::Type::VERTEX), ModelArray::Type::VERTEX },
    };
}

ParaGridIO::~ParaGridIO()
{
    for (auto& entry : openFiles) {
        entry.second.close();
    }
}

ModelState ParaGridIO::getModelState(const std::string& filePath)
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
    netCDF::NcGroup metaGroup(ncFile.getGroup(IStructure::metadataNodeName()));
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

    // Dimensions and DG components
    std::multimap<std::string, netCDF::NcDim> dimMap = dataGroup.getDims();
    for (auto entry : ModelArray::definedDimensions) {
        if (dimCompMap.count(entry.first) > 0)
            // TODO Assertions that DG in the file equals the compile time DG in the model. See #205
            continue;

        ModelArray::DimensionSpec& dimensionSpec = entry.second;
        netCDF::NcDim dim = dataGroup.getDim(dimensionSpec.name);
        if (entry.first == ModelArray::Dimension::Z) {
            // A special case, as the number of leves in the file might not be
            // the number that the selected ice thermodynamics requires.
            ModelArray::setDimension(entry.first, NZLevels::get());
        } else {
            ModelArray::setDimension(entry.first, dim.getSize());
        }
    }

    ModelState state;

    // Get all vars in the data group, and load them into a new ModelState

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
        ModelArray::Type newType = dimensionKeys.at(dimKey);
        state.data[varName] = ModelArray(newType);
        ModelArray& data = state.data.at(varName);
        data.resize();

        if (newType == ModelArray::Type::Z) {
            std::vector<size_t> startVector(ModelArray::nDimensions(newType), 0);
            std::vector<size_t> extentVector = ModelArray::dimensions(newType);
            std::reverse(extentVector.begin(), extentVector.end());
            var.getVar(startVector, extentVector, &data[0]);
        } else {
            var.getVar(&data[0]);
        }
    }
    ncFile.close();
    return state;
}

void ParaGridIO::dumpModelState(const ModelState& state, const ModelMetadata& metadata,
    const std::string& filePath, bool isRestart)
{

    if (isRestart) {
        netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);

        CommonRestartMetadata::writeStructureType(ncFile, metadata);
        netCDF::NcGroup metaGroup = ncFile.addGroup(IStructure::metadataNodeName());
        netCDF::NcGroup dataGroup = ncFile.addGroup(IStructure::dataNodeName());

        CommonRestartMetadata::writeRestartMetadata(metaGroup, metadata);

        // Dump the dimensions and number of components
        std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
        for (auto entry : ModelArray::definedDimensions) {
            ModelArray::Dimension dim = entry.first;
            size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim))
                                                   : dimSz = entry.second.length;
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
        // Reverse the order of the dimensions to translate between column-major ModelArray and
        // row-major netCDF
        for (auto& [type, v] : dimMap) {
            std::reverse(v.begin(), v.end());
        }
        // Everything that has components needs that dimension, too. This always varies fastest, and
        // so is last in the vector of dimensions.
        for (auto entry : dimCompMap) {
            dimMap.at(entry.second).push_back(ncFromMAMap.at(entry.first));
        }

        std::set<std::string> restartFields
            = { hiceName, ciceName, hsnowName, ticeName, maskName, coordsName }; // TODO and others
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
    } else {
        if (openFiles.count(filePath)) {
            // Append to an existing diagnostic output
            size_t nt = ++timeIndexByFile.at(filePath);
            netCDF::NcFile& ncFile = openFiles.at(filePath);

            netCDF::NcGroup metaGroup = ncFile.getGroup(IStructure::metadataNodeName());
            netCDF::NcGroup dataGroup = ncFile.getGroup(IStructure::dataNodeName());

            // Create references to the (existing) dimensions
            // Add the unlimited time dimension
            netCDF::NcDim timeDim = dataGroup.getDim(timeName);
            // All of the dimensions defined by the data at a particular timestep.
            std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
            for (auto entry : ModelArray::definedDimensions) {
                ModelArray::Dimension dim = entry.first;
                size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim))
                                                       : dimSz = entry.second.length;
                ncFromMAMap[dim] = dataGroup.getDim(entry.second.name);
            }

            // Also create the sets of dimensions to be connected to the data fields
            std::map<ModelArray::Type, std::vector<netCDF::NcDim>> dimMap;
            // Create the index and size arrays
            // The index arrays always start from zero, except in the first/time axis
            std::map<ModelArray::Type, std::vector<size_t>> indexArrays;
            std::map<ModelArray::Type, std::vector<size_t>> extentArrays;
            for (auto entry : ModelArray::typeDimensions) {
                ModelArray::Type type = entry.first;
                // No need to treat VERTEX arrays, they are assumed to be time constant
                if (type == ModelArray::Type::VERTEX)
                    continue;
                std::vector<netCDF::NcDim> ncDims;
                std::vector<size_t> indexArray;
                std::vector<size_t> extentArray;
                // Time dimension
                indexArray.push_back(nt);
                extentArray.push_back(1);
                // Other dimensions
                for (auto iter = entry.second.rbegin(); iter != entry.second.rend(); ++iter) {
                    ModelArray::Dimension& maDim = *iter;
                    ncDims.push_back(ncFromMAMap.at(maDim));
                    indexArray.push_back(0);
                    extentArray.push_back(ModelArray::definedDimensions.at(maDim).length);
                }
                dimMap[type] = ncDims;
                indexArrays[type] = indexArray;
                extentArrays[type] = extentArray;
            }
            // Everything that has components needs that dimension, too
            for (auto entry : dimCompMap) {
                if (entry.second == ModelArray::Type::VERTEX)
                    continue;
                dimMap.at(entry.second).push_back(ncFromMAMap.at(entry.first));
                indexArrays.at(entry.second).push_back(0);
                extentArrays.at(entry.second).push_back(ModelArray::nComponents(entry.second));
            }

            for (auto entry : state.data) {
                ModelArray::Type type = entry.second.getType();
                if (entry.first == maskName || type == ModelArray::Type::VERTEX)
                    continue;
                // Get the type, then relevant vector of NetCDF dimensions
                netCDF::NcVar var(dataGroup.getVar(entry.first));
                var.putVar(indexArrays.at(type), extentArrays.at(type), entry.second.getData());
            }
        } else {
            // Open a new diagnostic output file
            auto [it, success] = openFiles.try_emplace(filePath, filePath, netCDF::NcFile::replace);
            netCDF::NcFile& ncFile = openFiles.at(filePath);
            timeIndexByFile[filePath] = 0;

            CommonRestartMetadata::writeStructureType(ncFile, metadata);
            netCDF::NcGroup metaGroup = ncFile.addGroup(IStructure::metadataNodeName());
            netCDF::NcGroup dataGroup = ncFile.addGroup(IStructure::dataNodeName());

            CommonRestartMetadata::writeRestartMetadata(metaGroup, metadata);

            // Add the unlimited time dimension
            netCDF::NcDim timeDim = dataGroup.addDim(timeName);
            // All of the dimensions defined by the data at a particular timestep.
            std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
            for (auto entry : ModelArray::definedDimensions) {
                ModelArray::Dimension dim = entry.first;
                size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim))
                                                       : dimSz = entry.second.length;
                ncFromMAMap[dim] = dataGroup.addDim(entry.second.name, dimSz);
            }

            // Also create the sets of dimensions to be connected to the data fields
            std::map<ModelArray::Type, std::vector<netCDF::NcDim>> dimMap;
            // Create the index and size arrays
            // The index arrays always start from zero, except in the first/time axis
            std::map<ModelArray::Type, std::vector<size_t>> indexArrays;
            std::map<ModelArray::Type, std::vector<size_t>> extentArrays;
            for (auto entry : ModelArray::typeDimensions) {
                ModelArray::Type type = entry.first;
                std::vector<netCDF::NcDim> ncDims;
                // VERTEX arrays do not need time axes.
                std::vector<size_t> indexArray;
                std::vector<size_t> extentArray;
                if (type != ModelArray::Type::VERTEX) {
                    ncDims.push_back(timeDim);
                    indexArray.push_back(0);
                    extentArray.push_back(1UL);
                }
                for (ModelArray::Dimension& maDim : entry.second) {
                    ncDims.push_back(ncFromMAMap.at(maDim));
                    indexArray.push_back(0);
                    extentArray.push_back(ModelArray::definedDimensions.at(maDim).length);
                }
                dimMap[type] = ncDims;
                indexArrays[type] = indexArray;
                extentArrays[type] = extentArray;
            }
            // Everything that has components needs that dimension, too
            for (auto entry : dimCompMap) {
                dimMap.at(entry.second).push_back(ncFromMAMap.at(entry.first));
                indexArrays.at(entry.second).push_back(0);
                extentArrays.at(entry.second).push_back(ModelArray::nComponents(entry.second));
            }

            // Create a special timeless set of dimensions for the landmask
            std::vector<netCDF::NcDim> maskDims;
            for (ModelArray::Dimension& maDim :
                ModelArray::typeDimensions.at(ModelArray::Type::H)) {
                maskDims.push_back(ncFromMAMap.at(maDim));
            }
            std::vector<size_t> maskIndexes = { 0, 0 };
            std::vector<size_t> maskExtents
                = { ModelArray::definedDimensions
                          .at(ModelArray::typeDimensions.at(ModelArray::Type::H)[0])
                          .length,
                      ModelArray::definedDimensions
                          .at(ModelArray::typeDimensions.at(ModelArray::Type::H)[1])
                          .length };

            for (auto entry : state.data) {
                if (entry.first != maskName) {
                    // Get the type, then relevant vector of NetCDF dimensions
                    ModelArray::Type type = entry.second.getType();
                    std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
                    netCDF::NcVar var(dataGroup.addVar(entry.first, netCDF::ncDouble, ncDims));
                    var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
                    var.putVar(indexArrays.at(type), extentArrays.at(type), entry.second.getData());
                } else {
                    // Write the mask data
                    netCDF::NcVar var(dataGroup.addVar(maskName, netCDF::ncDouble, maskDims));
                    // No missing data
                    var.putVar(maskIndexes, maskExtents, entry.second.getData());
                }
            }
        }
    }
}

void ParaGridIO::close(const std::string& filePath)
{
    if (openFiles.count(filePath) > 0) {
        openFiles.at(filePath).close();
        openFiles.erase(openFiles.find(filePath));
        timeIndexByFile.erase(timeIndexByFile.find(filePath));
    }
}

} /* namespace Nextsim */
