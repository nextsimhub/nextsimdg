/*!
 * @file ParaGridIO.cpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ParaGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/FileCallbackCloser.hpp"
#include "include/Finalizer.hpp"
#include "include/MissingData.hpp"
#include "include/NZLevels.hpp"
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncException.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

#include <algorithm>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>

namespace Nextsim {

ParaGridIO::ParaGridIO(ParametricGrid& grid)
    : IParaGridIO(grid)
    , openFilesAndIndices(getOpenFilesAndIndices())
    , dimensionKeys({
          // clang-format off
          // Accept post-May 2024 (xdim, ydim, zdim) dimension names and pre-May 2024 (x, y, z)
        { "yx", ModelArray::Type::H },
        { "ydimxdim", ModelArray::Type::H },
        { "zyx", ModelArray::Type::Z },
        { "zdimydimxdim", ModelArray::Type::Z },
        { "yxdg_comp", ModelArray::Type::DG },
        { "ydimxdimdg_comp", ModelArray::Type::DG },
        { "yxdgstress_comp", ModelArray::Type::DGSTRESS },
        { "ydimxdimdgstress_comp", ModelArray::Type::DGSTRESS },
        { "ycgxcg", ModelArray::Type::CG },
        { "yvertexxvertexncoords", ModelArray::Type::VERTEX },
          // clang-format on
      })
    , isDG({
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
      })
    , dimCompMap({
          // clang-format off
        { ModelArray::componentMap.at(ModelArray::Type::DG), ModelArray::Type::DG },
        { ModelArray::componentMap.at(ModelArray::Type::DGSTRESS), ModelArray::Type::DGSTRESS },
        { ModelArray::componentMap.at(ModelArray::Type::VERTEX), ModelArray::Type::VERTEX },
          // clang-format on
      })
{
    if (doOnce()) {
        // Register the finalization function here
        Finalizer::atfinalUnique(closeAllFiles);
        // Since it should only ever run once, do further one-off initialization: allow distant
        // classes to close files via a callback.
        FileCallbackCloser::onClose(ParaGridIO::close);
        doOnce() = false;
    }
}

ParaGridIO::~ParaGridIO() = default;

#ifdef USE_MPI
ModelState ParaGridIO::getModelState(const std::string& filePath, ModelMetadata& metadata)
#else
ModelState ParaGridIO::getModelState(const std::string& filePath)
#endif
{
    ModelState state;

    try {
#ifdef USE_MPI
        netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::read, metadata.mpiComm);
#else
        netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
#endif
        netCDF::NcGroup metaGroup(ncFile.getGroup(IStructure::metadataNodeName()));
        netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

        // Dimensions and DG components
        std::multimap<std::string, netCDF::NcDim> dimMap = dataGroup.getDims();
        for (auto entry : ModelArray::definedDimensions) {
            auto dimType = entry.first;
            if (dimCompMap.count(dimType) > 0)
                // TODO Assertions that DG in the file equals the compile time DG in the model. See
                // #205
                continue;

            ModelArray::DimensionSpec& dimensionSpec = entry.second;
            // Find dimensions in the netCDF file by their name in the ModelArray details
            netCDF::NcDim dim = dataGroup.getDim(dimensionSpec.name);
            // Also check the old name
            if (dim.isNull()) {
                dim = dataGroup.getDim(dimensionSpec.altName);
            }
            // If we didn't find a dimension with the dimensions name or altName, throw.
            if (dim.isNull()) {
                throw std::out_of_range(
                    std::string("No netCDF dimension found corresponding to the dimension named ")
                    + dimensionSpec.name + std::string(" or ") + dimensionSpec.altName);
            }
            if (dimType == ModelArray::Dimension::Z) {
                // A special case, as the number of levels in the file might not be
                // the number that the selected ice thermodynamics requires.
#ifdef USE_MPI
                ModelArray::setDimension(dimType, NZLevels::get(), NZLevels::get(), 0);
#else
                ModelArray::setDimension(dimType, NZLevels::get());
#endif
            } else {
#ifdef USE_MPI
                auto dimName = dim.getName();
                size_t localLength = 0;
                size_t start = 0;
                if (dimType == ModelArray::Dimension::X) {
                    localLength = metadata.localExtentX;
                    start = metadata.localCornerX;
                } else if (dimType == ModelArray::Dimension::Y) {
                    localLength = metadata.localExtentY;
                    start = metadata.localCornerY;
                } else if (dimType == ModelArray::Dimension::XVERTEX) {
                    localLength = metadata.localExtentX + 1;
                    start = metadata.localCornerX;
                } else if (dimType == ModelArray::Dimension::YVERTEX) {
                    localLength = metadata.localExtentY + 1;
                    start = metadata.localCornerY;
                } else {
                    localLength = dim.getSize();
                    start = 0;
                }
                ModelArray::setDimension(dimType, dim.getSize(), localLength, start);
#else
                ModelArray::setDimension(dimType, dim.getSize());
#endif
            }
        }

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
                    std::string("No ModelArray::Type corresponds to the dimensional key ")
                    + dimKey);
            }
            ModelArray::Type type = dimensionKeys.at(dimKey);
            state.data[varName] = ModelArray(type);
            ModelArray& data = state.data.at(varName);
            data.resize();

            std::vector<size_t> start;
            std::vector<size_t> count;
            if (ModelArray::hasDoF(type)) {
                auto ncomps = data.nComponents();
                start.push_back(0);
                count.push_back(ncomps);
            }
            for (ModelArray::Dimension dt : ModelArray::typeDimensions.at(type)) {
                auto dim = ModelArray::definedDimensions.at(dt);
                start.push_back(dim.start);
                count.push_back(dim.localLength);
            }
            // dims are looped in [dg], x, y, [z] order so start and count
            // order must be reveresed to match order netcdf expects
            std::reverse(start.begin(), start.end());
            std::reverse(count.begin(), count.end());

            var.getVar(start, count, &data[0]);
        }
        ncFile.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filePath;
        throw std::runtime_error(ncWhat);
    }
    return state;
}

ModelState ParaGridIO::readForcingTimeStatic(
    const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath)
{
    ModelState state;

    try {
        netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
        netCDF::NcGroup metaGroup(ncFile.getGroup(IStructure::metadataNodeName()));
        netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

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
        const std::vector<ModelArray::Dimension>& dimensions
            = ModelArray::typeDimensions.at(ModelArray::Type::H);
        for (auto riter = dimensions.rbegin(); riter != dimensions.rend(); ++riter) {
            indexArray.push_back(0);
            extentArray.push_back(ModelArray::definedDimensions.at(*riter).localLength);
        }

        for (const std::string& varName : forcings) {
            netCDF::NcVar var = dataGroup.getVar(varName);
            state.data[varName] = ModelArray(ModelArray::Type::H);
            ModelArray& data = state.data.at(varName);
            data.resize();

            var.getVar(indexArray, extentArray, &data[0]);
        }
        ncFile.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filePath;
        throw std::runtime_error(ncWhat);
    }
    return state;
}

void ParaGridIO::dumpModelState(
    const ModelState& state, const ModelMetadata& metadata, const std::string& filePath)
{

#ifdef USE_MPI
    netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::replace, metadata.mpiComm);
#else
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);
#endif

    CommonRestartMetadata::writeStructureType(ncFile, metadata);
    netCDF::NcGroup metaGroup = ncFile.addGroup(IStructure::metadataNodeName());
    netCDF::NcGroup dataGroup = ncFile.addGroup(IStructure::dataNodeName());
    CommonRestartMetadata::writeRestartMetadata(metaGroup, metadata);

    // Dump the dimensions and number of components
    std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
    for (auto entry : ModelArray::definedDimensions) {
        ModelArray::Dimension dim = entry.first;
        size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim))
                                               : dimSz = entry.second.globalLength;
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
        sssName, maskName, coordsName, xName, yName, longitudeName, latitudeName, gridAzimuthName,
        uName, vName, damageName }; // TODO and others
    // If the above fields are found in the supplied ModelState, output them
    for (auto entry : state.data) {
        if (restartFields.count(entry.first)) {
            // Get the type, then relevant vector of NetCDF dimensions
            ModelArray::Type type = entry.second.getType();
            std::vector<size_t> start;
            std::vector<size_t> count;
            if (ModelArray::hasDoF(type)) {
                auto ncomps = entry.second.nComponents();
                start.push_back(0);
                count.push_back(ncomps);
            }
            for (ModelArray::Dimension dt : entry.second.typeDimensions.at(type)) {
                auto dim = entry.second.definedDimensions.at(dt);
                start.push_back(dim.start);
                count.push_back(dim.localLength);
            }
            // dims are looped in [dg], x, y, [z] order so start and count
            // order must be reveresed to match order netcdf expects
            std::reverse(start.begin(), start.end());
            std::reverse(count.begin(), count.end());

            std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
            netCDF::NcVar var(dataGroup.addVar(entry.first, netCDF::ncDouble, ncDims));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
            var.putVar(start, count, entry.second.getData());
        }
    }

    ncFile.close();
}

void ParaGridIO::writeDiagnosticTime(
    const ModelState& state, const ModelMetadata& meta, const std::string& filePath)
{
    bool isNew = openFilesAndIndices.count(filePath) <= 0;
    size_t nt = (isNew) ? 0 : ++openFilesAndIndices.at(filePath).second;
    if (isNew) {
        // Open a new file and emplace it in the map of open files.
        // Set the initial time to be zero (assigned above)
        // Piecewise construction is necessary to correctly construct the file handle/time index
        // pair
#ifdef USE_MPI
        openFilesAndIndices.emplace(std::piecewise_construct, std::make_tuple(filePath),
            std::forward_as_tuple(std::piecewise_construct,
                std::forward_as_tuple(filePath, netCDF::NcFile::replace, meta.mpiComm),
                std::forward_as_tuple(nt)));
#else
        openFilesAndIndices.emplace(std::piecewise_construct, std::make_tuple(filePath),
            std::forward_as_tuple(std::piecewise_construct,
                std::forward_as_tuple(filePath, netCDF::NcFile::replace),
                std::forward_as_tuple(nt)));
#endif
    }
    // Get the file handle
    NetCDFFileType& ncFile = openFilesAndIndices.at(filePath).first;

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
                                               : dimSz = entry.second.globalLength;
        ncFromMAMap[dim] = (isNew) ? dataGroup.addDim(entry.second.name, dimSz)
                                   : dataGroup.getDim(entry.second.name);
    }

    // Also create the sets of dimensions to be connected to the data fields
    std::map<ModelArray::Type, std::vector<netCDF::NcDim>> dimMap;
    // Create the index and size arrays
    std::map<ModelArray::Type, std::vector<size_t>> startMap;
    std::map<ModelArray::Type, std::vector<size_t>> countMap;
    for (auto entry : ModelArray::typeDimensions) {
        ModelArray::Type type = entry.first;
        std::vector<netCDF::NcDim> ncDims;
        std::vector<size_t> start;
        std::vector<size_t> count;

        // Everything that has components needs that dimension, too
        if (ModelArray::hasDoF(type)) {
            if (type == ModelArray::Type::VERTEX && !isNew)
                continue;
            auto ncomps = ModelArray::nComponents(type);
            auto dim = ModelArray::componentMap.at(type);
            ncDims.push_back(ncFromMAMap.at(dim));
            start.push_back(0);
            count.push_back(ncomps);
        }
        for (auto dt : entry.second) {
            auto dim = ModelArray::definedDimensions.at(dt);
            ncDims.push_back(ncFromMAMap.at(dt));
            start.push_back(dim.start);
            count.push_back(dim.localLength);
        }

        // Deal with VERTEX in each case
        // Add the time dimension for all types that are not VERTEX
        if (type != ModelArray::Type::VERTEX) {
            ncDims.push_back(timeDim);
            start.push_back(nt);
            count.push_back(1UL);
        } else if (!isNew) {
            // For VERTEX in an existing file, there is nothing more to be done
            continue;
        }

        std::reverse(ncDims.begin(), ncDims.end());
        std::reverse(start.begin(), start.end());
        std::reverse(count.begin(), count.end());

        dimMap[type] = ncDims;
        startMap[type] = start;
        countMap[type] = count;
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
                            .localLength,
            ModelArray::definedDimensions.at(ModelArray::typeDimensions.at(ModelArray::Type::H)[1])
                .localLength };
    }

    // Put the time axis variable
    std::vector<netCDF::NcDim> timeDimVec = { timeDim };
    netCDF::NcVar timeVar((isNew) ? dataGroup.addVar(timeName, netCDF::ncDouble, timeDimVec)
                                  : dataGroup.getVar(timeName));
    double secondsSinceEpoch = (meta.time() - TimePoint()).seconds();
#ifdef USE_MPI
    netCDF::setVariableCollective(timeVar, dataGroup);
#endif
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
#ifdef USE_MPI
            netCDF::setVariableCollective(var, dataGroup);
#endif
            var.putVar(maskIndexes, maskExtents, entry.second.getData());

        } else {
            std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
            // Get the variable object, either creating a new one or getting the existing one
            netCDF::NcVar var((isNew) ? dataGroup.addVar(entry.first, netCDF::ncDouble, ncDims)
                                      : dataGroup.getVar(entry.first));
            if (isNew)
                var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
#ifdef USE_MPI
            netCDF::setVariableCollective(var, dataGroup);
#endif
            var.putVar(startMap.at(type), countMap.at(type), entry.second.getData());
        }
    }
}

void ParaGridIO::close(const std::string& filePath)
{
    if (getOpenFilesAndIndices().count(filePath) > 0) {
        getOpenFilesAndIndices().at(filePath).first.close();
        getOpenFilesAndIndices().erase(filePath);
    }
}

void ParaGridIO::closeAllFiles()
{
    std::cout << "ParaGridIO::closeAllFiles: closing " << getOpenFilesAndIndices().size()
              << " files" << std::endl;
    while (getOpenFilesAndIndices().size() > 0) {
        close(getOpenFilesAndIndices().begin()->first);
    }
}

} /* namespace Nextsim */
