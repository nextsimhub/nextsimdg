/*!
 * @file ParaGridIO.cpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ParaGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/FileCallbackCloser.hpp"
#include "include/MissingData.hpp"
#include "include/NZLevels.hpp"
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

#include <algorithm>
#include <cstdlib>
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
std::map<std::string, netCDF::NcFile> ParaGridIO::openFiles;
std::map<std::string, size_t> ParaGridIO::timeIndexByFile;

void ParaGridIO::makeDimCompMap()
{
    dimCompMap = {
        { ModelArray::componentMap.at(ModelArray::Type::DG), ModelArray::Type::DG },
        { ModelArray::componentMap.at(ModelArray::Type::DGSTRESS), ModelArray::Type::DGSTRESS },
        { ModelArray::componentMap.at(ModelArray::Type::VERTEX), ModelArray::Type::VERTEX },
    };
    // Also initialize the static map of tables and register the atexit
    // function here, since it should only ever run once
    //    openFiles.clear();
    std::atexit(closeAllFiles);
    // Further one-off initialization: allow distant classes to close files via a callback.
    FileCallbackCloser::onClose(ParaGridIO::close);
}

ParaGridIO::~ParaGridIO() = default;

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
            // A special case, as the number of levels in the file might not be
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

ModelState ParaGridIO::readForcingTimeStatic(
    const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath)
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
    netCDF::NcGroup metaGroup(ncFile.getGroup(IStructure::metadataNodeName()));
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

    ModelState state;

    // Read the time axis
    netCDF::NcDim timeDim = dataGroup.getDim(timeName);
    // Read the time variable
    netCDF::NcVar timeVar = dataGroup.getVar(timeName);
    // Calculate the index of the largest time value on the axis below our target
    size_t targetTIndex;
    // Get the time axis as a vector
    std::vector<double> timeVec(timeDim.getSize());
    timeVar.getVar(timeVec.data());
    // Get the index of the largest TimePoint less than the target.
    targetTIndex = std::find_if(begin(timeVec), end(timeVec), [time](double t) {
        return (TimePoint() + Duration(t)) > time;
    }) - timeVec.begin();
    // Rather than the first that is greater than, get the last that is less
    // than or equal to. But don't go out of bounds.
    if (targetTIndex > 0)
        --targetTIndex;
    // ASSUME all forcings are HFields: finite volume fields on the same
    // grid as ice thickness
    std::vector<size_t> indexArray = { targetTIndex };
    std::vector<size_t> extentArray = { 1 };

    // Loop over the dimensions of H
    std::vector<ModelArray::Dimension>& dimensions
        = ModelArray::typeDimensions.at(ModelArray::Type::H);
    for (auto riter = dimensions.rbegin(); riter != dimensions.rend(); ++riter) {
        indexArray.push_back(0);
        extentArray.push_back(ModelArray::definedDimensions.at(*riter).length);
    }

    for (const std::string& varName : forcings) {
        netCDF::NcVar var = dataGroup.getVar(varName);
        state.data[varName] = ModelArray(ModelArray::Type::H);
        ModelArray& data = state.data.at(varName);
        data.resize();

        var.getVar(indexArray, extentArray, &data[0]);
    }
    ncFile.close();
    return state;
}

void ParaGridIO::dumpModelState(
    const ModelState& state, const ModelMetadata& metadata, const std::string& filePath)
{
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
        for (auto iter = entry.second.rbegin(); iter != entry.second.rend(); ++iter) {
            ModelArray::Dimension& maDim = *iter;
            ncDims.push_back(ncFromMAMap.at(maDim));
        }
        dimMap[type] = ncDims;
    }

    // Everything that has components needs that dimension, too. This always varies fastest, and so
    // is last in the vector of dimensions.
    for (auto entry : dimCompMap) {
        dimMap.at(entry.second).push_back(ncFromMAMap.at(entry.first));
    }

    std::set<std::string> restartFields = { hiceName, ciceName, hsnowName, ticeName, sstName,
        sssName, maskName, coordsName }; // TODO and others
    // Loop through either the above list (isRestart) or all provided fields(!isRestart)
    for (auto entry : state.data) {
        if (restartFields.count(entry.first)) {
            // Get the type, then relevant vector of NetCDF dimensions
            ModelArray::Type type = entry.second.getType();
            std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
            netCDF::NcVar var(dataGroup.addVar(entry.first, netCDF::ncDouble, ncDims));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value);
            var.putVar(entry.second.getData());
        }
    }

    ncFile.close();
}

void ParaGridIO::writeDiagnosticTime(
    const ModelState& state, const ModelMetadata& meta, const std::string& filePath)
{
    bool isNew = openFiles.count(filePath) <= 0;
    size_t nt = (isNew) ? 0 : ++timeIndexByFile.at(filePath);
    if (isNew) {
        // Open a new file and emplace it in the map of open files.
        openFiles.try_emplace(filePath, filePath, netCDF::NcFile::replace);
        // Set the initial time to be zero (assigned above)
        timeIndexByFile[filePath] = nt;
    }
    // Get the file handle
    netCDF::NcFile& ncFile = openFiles.at(filePath);

    // Get the netCDF groups, creating them if necessary
    netCDF::NcGroup metaGroup = (isNew) ? ncFile.addGroup(IStructure::metadataNodeName())
                                        : ncFile.getGroup(IStructure::metadataNodeName());
    netCDF::NcGroup dataGroup = (isNew) ? ncFile.addGroup(IStructure::dataNodeName())
                                        : ncFile.getGroup(IStructure::dataNodeName());

    if (isNew) {
        // Write the common structure and time metadata
        CommonRestartMetadata::writeStructureType(ncFile, meta);
        CommonRestartMetadata::writeRestartMetadata(metaGroup, meta);
    }
    // Get the unlimited time dimension, creating it if necessary
    netCDF::NcDim timeDim = (isNew) ? dataGroup.addDim(timeName) : dataGroup.getDim(timeName);

    // All of the dimensions defined by the data at a particular timestep.
    std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
    for (auto entry : ModelArray::definedDimensions) {
        ModelArray::Dimension dim = entry.first;
        size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim))
                                               : dimSz = entry.second.length;
        ncFromMAMap[dim] = (isNew) ? dataGroup.addDim(entry.second.name, dimSz)
                                   : dataGroup.getDim(entry.second.name);
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
        std::vector<size_t> indexArray;
        std::vector<size_t> extentArray;

        // Deal with VERTEX in each case
        // Add the time dimension for all types that are not VERTEX
        if (type != ModelArray::Type::VERTEX) {
            ncDims.push_back(timeDim);
            indexArray.push_back(nt);
            extentArray.push_back(1UL);
        } else if (!isNew) {
            // For VERTEX in an existing file, there is nothing more to be done
            continue;
        }
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
        // Skip VERTEX fields on existing files
        if (entry.second == ModelArray::Type::VERTEX && !isNew)
            continue;
        dimMap.at(entry.second).push_back(ncFromMAMap.at(entry.first));
        indexArrays.at(entry.second).push_back(0);
        extentArrays.at(entry.second).push_back(ModelArray::nComponents(entry.second));
    }

    // Create a special timeless set of dimensions for the landmask
    std::vector<netCDF::NcDim> maskDims;
    std::vector<size_t> maskIndexes;
    std::vector<size_t> maskExtents;
    if (isNew) {
        for (ModelArray::Dimension& maDim : ModelArray::typeDimensions.at(ModelArray::Type::H)) {
            maskDims.push_back(ncFromMAMap.at(maDim));
        }
        maskIndexes = { 0, 0 };
        maskExtents = { ModelArray::definedDimensions
                            .at(ModelArray::typeDimensions.at(ModelArray::Type::H)[0])
                            .length,
            ModelArray::definedDimensions.at(ModelArray::typeDimensions.at(ModelArray::Type::H)[1])
                .length };
    }

    // Put the time axis variable
    std::vector<netCDF::NcDim> timeDimVec = { timeDim };
    netCDF::NcVar timeVar((isNew) ? dataGroup.addVar(timeName, netCDF::ncDouble, timeDimVec)
                                  : dataGroup.getVar(timeName));
    double secondsSinceEpoch = (meta.time() - TimePoint()).seconds();
    timeVar.putVar({ nt }, { 1 }, &secondsSinceEpoch);

    // Write the data
    for (auto entry : state.data) {
        ModelArray::Type type = entry.second.getType();
        // Skip timeless fields (mask, coordinates) on existing files
        if (!isNew && (entry.first == maskName || type == ModelArray::Type::VERTEX))
            continue;
        if (entry.first == maskName) {
            // Land mask in a new file (since it was skipped above in existing files)
            netCDF::NcVar var(dataGroup.addVar(maskName, netCDF::ncDouble, maskDims));
            // No missing data
            var.putVar(maskIndexes, maskExtents, entry.second.getData());

        } else {
            std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
            // Get the variable object, either creating a new one or getting the existing one
            netCDF::NcVar var((isNew) ? dataGroup.addVar(entry.first, netCDF::ncDouble, ncDims)
                                      : dataGroup.getVar(entry.first));
            if (isNew)
                var.putAtt(mdiName, netCDF::ncDouble, MissingData::value);

            var.putVar(indexArrays.at(type), extentArrays.at(type), entry.second.getData());
        }
    }
}

void ParaGridIO::close(const std::string& filePath)
{
    if (openFiles.count(filePath) > 0) {
        openFiles.at(filePath).close();
        openFiles.erase(filePath);
        timeIndexByFile.erase(filePath);
    }
}

void ParaGridIO::closeAllFiles()
{
    size_t closedFiles = 0;
    for (const auto& [name, handle] : openFiles) {
        if (!handle.isNull()) {
            close(name);
            closedFiles++;
        }
        /* If the following break is not checked for and performed, for some
         * reason the iteration will continue to iterate over invalid
         * string/NcFile pairs. */
        if (closedFiles >= openFiles.size())
            break;
    }
}

} /* namespace Nextsim */
