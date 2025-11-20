/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ParaGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/FileCallbackCloser.hpp"
#include "include/Finalizer.hpp"
#ifdef USE_MPI
#include "include/Halo.hpp"
#include "include/ModelMPI.hpp"
#endif
#include "include/MissingData.hpp"
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncException.h>
#include <ncFile.h>
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
        { "yxdg_comp", ModelArray::Type::DG },
        { "ydimxdimdg_comp", ModelArray::Type::DG },
        { "yxdgstress_comp", ModelArray::Type::DGSTRESS },
        { "ydimxdimdgstress_comp", ModelArray::Type::DGSTRESS },
        { "y_cgx_cg", ModelArray::Type::CG },
        { "yvertexxvertexncoords", ModelArray::Type::VERTEX },
          // clang-format on
      })
    , isDG({
          // clang-format off
        { ModelArray::Dimension::X, false },
        { ModelArray::Dimension::Y, false },
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
    static bool doneOnce = doOnce();
}

bool ParaGridIO::doOnce()
{
    // Register the finalization function here
    Finalizer::registerUnique(closeAllFiles);
    // Since it should only ever run once, do further one-off initialization: allow distant
    // classes to close files via a callback.
    FileCallbackCloser::onClose(close);

    return true;
}

ParaGridIO::~ParaGridIO() = default;

ModelState ParaGridIO::getModelState(const std::string& filePath)
{
    ModelState state;

    try {
#ifdef USE_MPI
        auto& modelMPI = ModelMPI::getInstance();
        netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::read, modelMPI.getComm());
#else
        netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
#endif

        // Dimensions and DG components
        std::multimap<std::string, netCDF::NcDim> dimMap = ncFile.getDims();
        for (auto entry : ModelArray::definedDimensions) {
            auto dimType = entry.first;
            if (dimCompMap.count(dimType) > 0)
                // TODO Assertions that DG in the file equals the compile time DG in the model. See
                // #205
                continue;

            ModelArray::DimensionSpec& dimensionSpec = entry.second;
            // Find dimensions in the netCDF file by their name in the ModelArray details
            netCDF::NcDim dim = ncFile.getDim(dimensionSpec.name);
            // Also check the old name
            if (dim.isNull()) {
                dim = ncFile.getDim(dimensionSpec.altName);
            }
            // If we didn't find a dimension with the dimensions name or altName, throw.
            if (dim.isNull()) {
                throw std::out_of_range(
                    std::string("No netCDF dimension found corresponding to the dimension named ")
                    + dimensionSpec.name + std::string(" or ") + dimensionSpec.altName);
            }
#ifdef USE_MPI
            auto dimName = dim.getName();
            size_t localLength = 0;
            size_t start = 0;
            auto& metadata = ModelMetadata::getInstance();
            if (dimType == ModelArray::Dimension::X) {
                localLength = metadata.getLocalExtentX();
                start = metadata.getLocalCornerX();
            } else if (dimType == ModelArray::Dimension::Y) {
                localLength = metadata.getLocalExtentY();
                start = metadata.getLocalCornerY();
            } else if (dimType == ModelArray::Dimension::XVERTEX) {
                localLength = metadata.getLocalExtentX() + 1;
                start = metadata.getLocalCornerX();
            } else if (dimType == ModelArray::Dimension::YVERTEX) {
                localLength = metadata.getLocalExtentY() + 1;
                start = metadata.getLocalCornerY();
            } else {
                localLength = dim.getSize();
                start = 0;
            }
            // globalLength doesnt need to be padded with halo cells but localLength does
            // setDimension(dim, globalLength, localLength, start)
            ModelArray::setDimension(
                dimType, dim.getSize(), localLength + 2 * Halo::haloWidth, start);
#else
            ModelArray::setDimension(dimType, dim.getSize());
#endif
        }

        // Get all valid variables and load them into a new ModelState

        for (auto entry : ncFile.getVars()) {
            const std::string& varName = entry.first;
            netCDF::NcVar& var = entry.second;
            // Determine the type from the dimensions
            std::vector<netCDF::NcDim> varDims = var.getDims();
            std::string dimKey = "";
            for (netCDF::NcDim& dim : varDims) {
                dimKey += dim.getName();
            }
            // Skip invalid dimension keys
            if (!dimensionKeys.count(dimKey)) {
                continue;
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

            using Dim = ModelArray::Dimension;
            for (Dim dt : ModelArray::typeDimensions.at(type)) {
                auto dim = ModelArray::definedDimensions.at(dt);
                size_t startIndex = dim.start;
                size_t localLength = dim.localLength;
#ifdef USE_MPI
                // Halo cells (which only exist in the lateral direction) are not included in netCDF
                // files. Only read/write the inner block
                if (Halo::isDimLateral(dt)) {
                    localLength = localLength - 2 * Halo::haloWidth;
                }
#endif
                start.push_back(startIndex);
                count.push_back(localLength);
            }
            // dims are looped in [dg], x, y, [z] order so start and count
            // order must be reversed to match order netcdf expects
            std::reverse(start.begin(), start.end());
            std::reverse(count.begin(), count.end());

#ifdef USE_MPI
            Halo halo(data);
            // create and allocate temporary Eigen array
            ModelArray::DataType tempData;
            tempData.resize(halo.getInnerSize(), data.nComponents());
            // populate temp Eigen array with data from netCDF file
            var.getVar(start, count, tempData.data());
            // populate inner block of modelarray with data from tempData
            halo.setInnerBlock(tempData, data.getDataRef());
            halo.exchangeHalos(data.getDataRef());
#else
            var.getVar(start, count, &data[0]);
#endif
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

        // Read the time axis
        netCDF::NcDim timeDim = ncFile.getDim(timeName);
        // Read the time variable
        netCDF::NcVar timeVar = ncFile.getVar(timeName);
        // Calculate the index of the largest time value on the axis below our target
        size_t targetTIndex;
        // Get the time axis as a vector
        std::vector<double> timeVec(timeDim.getSize());
        timeVar.getVar(timeVec.data());
        // Get the index of the largest TimePoint less than the target.
        targetTIndex = std::find_if(std::begin(timeVec), std::end(timeVec), [time](double t) {
            return (TimePoint() + Duration(t)) > time;
        }) - timeVec.begin();
        // Rather than the first that is greater than, get the last that is less
        // than or equal to. But don't go out of bounds.
        if (targetTIndex > 0)
            --targetTIndex;
        // ASSUME all forcings are HFields: finite volume fields on the same
        // grid as ice thickness
        std::vector<size_t> indexArray;
        std::vector<size_t> extentArray;

        // Loop over the dimensions of H
        using Dim = ModelArray::Dimension;
        for (Dim dt : ModelArray::typeDimensions.at(ModelArray::Type::H)) {
            auto dim = ModelArray::definedDimensions.at(dt);
            auto startIndex = dim.start;
            auto localLength = dim.localLength;
#ifdef USE_MPI
            // Halo cells (which only exist in the lateral direction) are not included in netCDF
            // files. Only read/write the inner block
            if (Halo::isDimLateral(dt)) {
                localLength = localLength - 2 * Halo::haloWidth;
            }
#endif
            indexArray.push_back(startIndex);
            extentArray.push_back(localLength);
        }

        indexArray.push_back(targetTIndex);
        extentArray.push_back(1);
        std::reverse(indexArray.begin(), indexArray.end());
        std::reverse(extentArray.begin(), extentArray.end());

        auto availableForcings = ncFile.getVars();
        for (const std::string& varName : forcings) {
            // Don't try to read non-existent data
            if (!availableForcings.count(varName)) {
                continue;
            }
            netCDF::NcVar var = ncFile.getVar(varName);
            state.data[varName] = ModelArray(ModelArray::Type::H);
            ModelArray& data = state.data.at(varName);
            data.resize();

#ifdef USE_MPI
            Halo halo(data);
            // create and allocate temporary Eigen array
            ModelArray::DataType tempData;
            tempData.resize(halo.getInnerSize(), data.nComponents());
            // populate temp Eigen array with data from netCDF file
            var.getVar(indexArray, extentArray, tempData.data());
            // populate inner block of modelarray with data from tempData
            halo.setInnerBlock(tempData, data.getDataRef());
            halo.exchangeHalos(data.getDataRef());
#else
            var.getVar(indexArray, extentArray, &data[0]);
#endif
        }
        ncFile.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filePath;
        throw std::runtime_error(ncWhat);
    }
    return state;
}

void ParaGridIO::dumpModelState(const ModelState& state, const std::string& filePath)
{

#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance();
    netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::replace, modelMPI.getComm());
#else
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);
#endif

    CommonRestartMetadata::writeStructureType(ncFile);
    CommonRestartMetadata::writeRestartMetadata(ncFile);

    // Dump the dimensions and number of components
    std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
    for (auto entry : ModelArray::definedDimensions) {
        ModelArray::Dimension dim = entry.first;
        size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim))
                                               : dimSz = entry.second.localLength;
        ncFromMAMap[dim] = ncFile.addDim(entry.second.name, dimSz);
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

    // Assume that all fields in the supplied ModelState are necessary, and so write them to
    // file.
    for (auto entry : state.data) {
        // Get the type, then relevant vector of NetCDF dimensions
        ModelArray::Type type = entry.second.getType();
        std::vector<size_t> start;
        std::vector<size_t> count;
        if (ModelArray::hasDoF(type)) {
            auto ncomps = entry.second.nComponents();
            start.push_back(0);
            count.push_back(ncomps);
        }
        using Dim = ModelArray::Dimension;
        for (ModelArray::Dimension dt : entry.second.typeDimensions.at(type)) {
            auto dim = entry.second.definedDimensions.at(dt);
            size_t localLength = dim.localLength;
#ifdef USE_MPI
            // Halo cells (which only exist in the lateral direction) are not included in netCDF
            // files. Only read/write the inner block
            if (Halo::isDimLateral(dt)) {
                localLength = localLength - 2 * Halo::haloWidth;
            }
#endif
            start.push_back(dim.start);
            count.push_back(localLength);
        }
        // dims are looped in [dg], x, y, [z] order so start and count
        // order must be reveresed to match order netcdf expects
        std::reverse(start.begin(), start.end());
        std::reverse(count.begin(), count.end());

        std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
        netCDF::NcVar var(ncFile.addVar(entry.first, netCDF::ncDouble, ncDims));
        var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());

#ifdef USE_MPI
        auto& data = entry.second;
        Halo halo(data);
        ModelArray::DataType tempData;
        tempData.resize(halo.getInnerSize(), data.nComponents());
        halo.getInnerBlock(data.getDataRef(), tempData);
        var.putVar(start, count, tempData.data());
#else
        var.putVar(start, count, entry.second.getData());
#endif
    }

    ncFile.close();
}

void ParaGridIO::writeDiagnosticTime(const ModelState& state, const std::string& filePath)
{
    bool isNew = openFilesAndIndices.count(filePath) <= 0;
    size_t nt = (isNew) ? 0 : ++openFilesAndIndices.at(filePath).second;
    if (isNew) {
        // Open a new file and emplace it in the map of open files.
        // Set the initial time to be zero (assigned above)
        // Piecewise construction is necessary to correctly construct the file handle/time index
        // pair
#ifdef USE_MPI
        auto& modelMPI = ModelMPI::getInstance();
        openFilesAndIndices.emplace(std::piecewise_construct, std::make_tuple(filePath),
            std::forward_as_tuple(std::piecewise_construct,
                std::forward_as_tuple(filePath, netCDF::NcFile::replace, modelMPI.getComm()),
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

    if (isNew) {
        // Write the common structure and time metadata
        CommonRestartMetadata::writeStructureType(ncFile);
    }
    // Get the unlimited time dimension, creating it if necessary
    netCDF::NcDim timeDim = (isNew) ? ncFile.addDim(timeName) : ncFile.getDim(timeName);

    // All of the dimensions defined by the data at a particular timestep.
    std::map<ModelArray::Dimension, netCDF::NcDim> ncFromMAMap;
    for (auto entry : ModelArray::definedDimensions) {
        ModelArray::Dimension dim = entry.first;
        size_t dimSz = (dimCompMap.count(dim)) ? ModelArray::nComponents(dimCompMap.at(dim))
                                               : dimSz = entry.second.globalLength;
        ncFromMAMap[dim]
            = (isNew) ? ncFile.addDim(entry.second.name, dimSz) : ncFile.getDim(entry.second.name);
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
        using Dim = ModelArray::Dimension;
        for (auto dt : entry.second) {
            auto dim = ModelArray::definedDimensions.at(dt);
            auto localLength = dim.localLength;
#ifdef USE_MPI
            // Halo cells (which only exist in the lateral direction) are not included in netCDF
            // files. Only read/write the inner block
            if (Halo::isDimLateral(dt)) {
                localLength = localLength - 2 * Halo::haloWidth;
            }
#endif
            ncDims.push_back(ncFromMAMap.at(dt));
            start.push_back(dim.start);
            count.push_back(localLength);
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
        for (auto dt : ModelArray::typeDimensions.at(ModelArray::Type::H)) {
            auto dim = ModelArray::definedDimensions.at(dt);
            auto localLength = dim.localLength;
#ifdef USE_MPI
            // Halo cells (which only exist in the lateral direction) are not included in netCDF
            // files. Only read/write the inner block
            if (Halo::isDimLateral(dt)) {
                localLength = localLength - 2 * Halo::haloWidth;
            }
#endif
            maskIndexes.push_back(0);
            maskExtents.push_back(localLength);
        }
    }

    // Put the time axis variable
    std::vector<netCDF::NcDim> timeDimVec = { timeDim };
    netCDF::NcVar timeVar(
        (isNew) ? ncFile.addVar(timeName, netCDF::ncDouble, timeDimVec) : ncFile.getVar(timeName));
    auto& metadata = ModelMetadata::getInstance();
    double secondsSinceEpoch = (metadata.time() - TimePoint()).seconds();
#ifdef USE_MPI
    netCDF::setVariableCollective(timeVar, ncFile);
#endif
    timeVar.putVar({ nt }, { 1 }, &secondsSinceEpoch);

    if (isNew)
        timeVar.putAtt("units", "seconds since 1970-01-01 00:00:00");

    // Write the data
    for (auto entry : state.data) {
        ModelArray::Type type = entry.second.getType();
        // Skip timeless fields (mask, coordinates) on existing files
        if (!isNew && (entry.first == maskName || type == ModelArray::Type::VERTEX))
            continue;
        if (entry.first == maskName) {
            // Land mask in a new file (since it was skipped above in existing files)
            netCDF::NcVar var(ncFile.addVar(maskName, netCDF::ncDouble, maskDims));
            // No missing data
#ifdef USE_MPI
            netCDF::setVariableCollective(var, ncFile);
            auto& data = entry.second;
            Halo halo(data);
            ModelArray::DataType tempData;
            tempData.resize(halo.getInnerSize(), data.nComponents());
            halo.getInnerBlock(data.getDataRef(), tempData);
            var.putVar(maskIndexes, maskExtents, tempData.data());
#else
            var.putVar(maskIndexes, maskExtents, entry.second.getData());
#endif

        } else {
            std::vector<netCDF::NcDim>& ncDims = dimMap.at(type);
            // Get the variable object, either creating a new one or getting the existing one
            netCDF::NcVar var((isNew) ? ncFile.addVar(entry.first, netCDF::ncDouble, ncDims)
                                      : ncFile.getVar(entry.first));
            if (isNew)
                var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
#ifdef USE_MPI
            netCDF::setVariableCollective(var, ncFile);
            auto& data = entry.second;
            Halo halo(data);
            ModelArray::DataType tempData;
            tempData.resize(halo.getInnerSize(), data.nComponents());
            halo.getInnerBlock(data.getDataRef(), tempData);
            var.putVar(startMap.at(type), countMap.at(type), tempData.data());
#else
            var.putVar(startMap.at(type), countMap.at(type), entry.second.getData());
#endif
        }
    }

    // Flush buffer to disc for monitoring and diagnostic purposes.
    ncFile.sync();
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
