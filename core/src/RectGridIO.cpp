/*!
 * @file RectGridIO.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#include "include/RectGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/IStructure.hpp"
#include "include/MissingData.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelState.hpp"
#include "include/NZLevels.hpp"
#include "include/gridNames.hpp"

#ifdef USE_MPI
#include "include/ParallelNetcdfFile.hpp"
#endif

#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncVar.h>

#include <algorithm>
#include <list>
#include <map>
#include <string>
#include <vector>

namespace Nextsim {

#ifdef USE_MPI
void dimensionSetter(const netCDF::NcGroup& dataGroup, const std::string& fieldName,
    ModelArray::Type type, ModelMetadata& metadata)
{
    size_t nDims = dataGroup.getVar(fieldName).getDimCount();
    ModelArray::MultiDim dims;
    dims.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        dims[d] = dataGroup.getVar(fieldName).getDim(d).getSize();
    }
    // The dimensions in the netCDF are in the reverse order compared to ModelArray
    std::reverse(dims.begin(), dims.end());
    // A special case for Type::Z: use NZLevels for the third dimension
    if (type == ModelArray::Type::Z)
        dims[2] = NZLevels::get();
    // Replace X, Y dimensions with local extends
    dims[0] = metadata.localExtentX;
    dims[1] = metadata.localExtentY;
    ModelArray::setDimensions(type, dims);
}
#else
void dimensionSetter(
    const netCDF::NcGroup& dataGroup, const std::string& fieldName, ModelArray::Type type)
{
    size_t nDims = dataGroup.getVar(fieldName).getDimCount();
    ModelArray::MultiDim dims;
    dims.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        dims[d] = dataGroup.getVar(fieldName).getDim(d).getSize();
    }
    // The dimensions in the netCDF are in the reverse order compared to ModelArray
    std::reverse(dims.begin(), dims.end());
    // A special case for Type::Z: use NZLevels for the third dimension
    if (type == ModelArray::Type::Z)
        dims[2] = NZLevels::get();
    ModelArray::setDimensions(type, dims);
}
#endif

#ifdef USE_MPI
ModelState RectGridIO::getModelState(const std::string& filePath, ModelMetadata& metadata)
#else
ModelState RectGridIO::getModelState(const std::string& filePath)
#endif
{
    ModelState state;
#ifdef USE_MPI
    netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::read, metadata.mpiComm);
#else
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
#endif
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

#ifdef USE_MPI
    // Get the sizes of the four types of field
    // HField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::H, metadata);
    // UField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::U, metadata);
    // VField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::V, metadata);
    // ZField from tice
    dimensionSetter(dataGroup, ticeName, ModelArray::Type::Z, metadata);
#else
    // Get the sizes of the four types of field
    // HField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::H);
    // UField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::U);
    // VField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::V);
    // ZField from tice
    dimensionSetter(dataGroup, ticeName, ModelArray::Type::Z);
#endif

#ifdef USE_MPI
    // Set the origins and extensions for reading 2D data based
    // on MPI decomposition
    std::vector<size_t> start(2);
    std::vector<size_t> size(2);
    start[0] = metadata.localCornerY;
    start[1] = metadata.localCornerX;
    size[0] = metadata.localExtentY;
    size[1] = metadata.localExtentX;
#else
    std::vector<size_t> start = { 0, 0 };
    std::vector<size_t> size = ModelArray::dimensions(ModelArray::Type::H);
    std::reverse(size.begin(), size.end());
#endif
    state.data[maskName] = ModelArray::HField();
    dataGroup.getVar(maskName).getVar(start, size, &state.data[maskName][0]);
    state.data[hiceName] = ModelArray::HField();
    dataGroup.getVar(hiceName).getVar(start, size, &state.data[hiceName][0]);
    state.data[ciceName] = ModelArray::HField();
    dataGroup.getVar(ciceName).getVar(start, size, &state.data[ciceName][0]);
    state.data[hsnowName] = ModelArray::HField();
    dataGroup.getVar(hsnowName).getVar(start, size, &state.data[hsnowName][0]);
    // coordinates on the H grid
    if (dataGroup.getVars().count(xName) > 0) {
        state.data[xName] = ModelArray::HField();
        dataGroup.getVar(xName).getVar(start, size, &state.data[xName][0]);
        state.data[yName] = ModelArray::HField();
        dataGroup.getVar(yName).getVar(start, size, &state.data[yName][0]);
    } else {
        state.data[longitudeName] = ModelArray::HField();
        dataGroup.getVar(longitudeName).getVar(start, size, &state.data[longitudeName][0]);
        state.data[latitudeName] = ModelArray::HField();
        dataGroup.getVar(latitudeName).getVar(start, size, &state.data[latitudeName][0]);
    }

    // Z direction is outside MPI ifdef as the domain is never decomposed in this direction

    // Since the ZFierld might not have the same dimensions as the tice field
    // in the file, a little more work is required.
    state.data[ticeName] = ModelArray::ZField();
    start.insert(start.begin(), 0);
    std::vector<size_t> zArrayDims = ModelArray::dimensions(ModelArray::Type::Z);
    std::reverse(zArrayDims.begin(), zArrayDims.end());
    dataGroup.getVar(ticeName).getVar(start, zArrayDims, &state.data[ticeName][0]);

    ncFile.close();
    return state;
}

void RectGridIO::dumpModelState(const ModelState& state, const ModelMetadata& metadata,
    const std::string& filePath, bool isRestart) const
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
    typedef ModelArray::Type Type;

    int nx = ModelArray::dimensions(Type::H)[0];
    int ny = ModelArray::dimensions(Type::H)[1];
    int nz = ModelArray::dimensions(Type::Z)[2];

    std::vector<std::string> dimensionNames
        = { "xdim", "ydim", "zdim", "t", "component", "u", "v", "w" };

    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
#ifdef USE_MPI
    netCDF::NcDim xDim = dataGroup.addDim(dimensionNames[0], metadata.globalExtentX);
    netCDF::NcDim yDim = dataGroup.addDim(dimensionNames[1], metadata.globalExtentY);
#else
    netCDF::NcDim xDim = dataGroup.addDim(dimensionNames[0], nx);
    netCDF::NcDim yDim = dataGroup.addDim(dimensionNames[1], ny);
#endif
    netCDF::NcDim zDim = dataGroup.addDim(dimensionNames[2], nz);
    std::vector<netCDF::NcDim> dims2 = { yDim, xDim };
    std::vector<netCDF::NcDim> dims3 = { zDim, yDim, xDim };
#ifdef USE_MPI
    // Set the origins and extensions for reading 3D data based
    // on MPI decomposition
    std::vector<size_t> start3 = { 0, static_cast<size_t>(metadata.localCornerY),
        static_cast<size_t>(metadata.localCornerX) };
    std::vector<size_t> size3 = { static_cast<size_t>(nz),
        static_cast<size_t>(metadata.localExtentY), static_cast<size_t>(metadata.localExtentX) };
    // Set the origins and extensions for reading 2D data based
    // on MPI decomposition
    std::vector<size_t> start2(start3.begin() + 1, start3.end());
    std::vector<size_t> size2(size3.begin() + 1, size3.end());
#endif

    for (const auto entry : state.data) {
        const std::string& name = entry.first;
        if (entry.second.getType() == ModelArray::Type::H && entry.second.trueSize() > 0) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims2));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value);
#ifdef USE_MPI
            var.putVar(start2, size2, entry.second.getData());
#else
            var.putVar(entry.second.getData());
#endif
        } else if (entry.second.getType() == ModelArray::Type::Z && entry.second.trueSize() > 0) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims3));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value);
#ifdef USE_MPI
            var.putVar(start3, size3, entry.second.getData());
#else
            var.putVar(entry.second.getData());
#endif
        }
    }

    ncFile.close();
}

} /* namespace Nextsim */
