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
#include "include/RectangularGrid.hpp"
#include "include/gridNames.hpp"

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
    dims[0] = metadata.localExtentX;
    dims[1] = metadata.localExtentY;
    for (size_t d = 2; d < nDims; ++d) {
        dims[d] = dataGroup.getVar(fieldName).getDim(d).getSize();
    }
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
ModelState RectGridIO::getModelState(
    const std::string& filePath, const std::string& partitionFile, ModelMetadata& metadata)
#else
ModelState RectGridIO::getModelState(const std::string& filePath)
#endif
{
    ModelState state;
#ifdef USE_MPI
    readPartitionData(partitionFile, metadata);
    netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::read, metadata.mpiComm);
#else
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
#endif
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

#ifdef USE_MPI
    // Get the sizes of the four types of field
    // HField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::H, metadata);
    // UField from hice        TODO replace with u velocity once it is present
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::U, metadata);
    // VField from hice        TODO replace with v velocity once it is present
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::V, metadata);
    // ZField from tice
    dimensionSetter(dataGroup, ticeName, ModelArray::Type::Z, metadata);

    // Set the origins and extensions for reading 2D data based
    // on MPI decomposition
    std::vector<size_t> start(2);
    std::vector<size_t> size(2);
    start[0] = metadata.localCornerX;
    start[1] = metadata.localCornerY;
    size[0] = metadata.localExtentX;
    size[1] = metadata.localExtentY;

    state.data[maskName] = ModelArray::HField();
    dataGroup.getVar(maskName).getVar(start, size, &state.data[maskName][0]);
    state.data[hiceName] = ModelArray::HField();
    dataGroup.getVar(hiceName).getVar(start, size, &state.data[hiceName][0]);
    state.data[ciceName] = ModelArray::HField();
    dataGroup.getVar(ciceName).getVar(start, size, &state.data[ciceName][0]);
    state.data[hsnowName] = ModelArray::HField();
    dataGroup.getVar(hsnowName).getVar(start, size, &state.data[hsnowName][0]);

    // The domain is not decomposed in z direction so set the extend in this direction
    // to full range
    start.push_back(0);
    size.push_back(dataGroup.getDim("nLayers").getSize());

    state.data[ticeName] = ModelArray::ZField();
    dataGroup.getVar(ticeName).getVar(start, size, &state.data[ticeName][0]);

#else
    // Get the sizes of the four types of field
    // HField from hice
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::H);
    // UField from hice        TODO replace with u velocity once it is present
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::U);
    // VField from hice        TODO replace with v velocity once it is present
    dimensionSetter(dataGroup, hiceName, ModelArray::Type::V);
    // ZField from tice
    dimensionSetter(dataGroup, ticeName, ModelArray::Type::Z);

    state.data[maskName] = ModelArray::HField();
    dataGroup.getVar(maskName).getVar(&state.data[maskName][0]);
    state.data[hiceName] = ModelArray::HField();
    dataGroup.getVar(hiceName).getVar(&state.data[hiceName][0]);
    state.data[ciceName] = ModelArray::HField();
    dataGroup.getVar(ciceName).getVar(&state.data[ciceName][0]);
    state.data[hsnowName] = ModelArray::HField();
    dataGroup.getVar(hsnowName).getVar(&state.data[hsnowName][0]);
    // Since the ZFierld might not have the same dimensions as the tice field
    // in the file, a little more work is required.
    state.data[ticeName] = ModelArray::ZField();
    std::vector<size_t> startVector = { 0, 0, 0 };
    std::vector<size_t> zArrayDims = ModelArray::dimensions(ModelArray::Type::Z);
    std::reverse(zArrayDims.begin(), zArrayDims.end());
    dataGroup.getVar(ticeName).getVar(startVector, zArrayDims, &state.data[ticeName][0]);
#endif

    ncFile.close();
    return state;
}

#ifdef USE_MPI
void RectGridIO::readPartitionData(const std::string& partitionFile, ModelMetadata& metadata)
{
    static const std::string bboxName = "bounding_boxes";

    netCDF::NcFile ncFile(partitionFile, netCDF::NcFile::read);
    int sizes = ncFile.getDim("L").getSize();
    int nBoxes = ncFile.getDim("P").getSize();
    if (nBoxes != metadata.mpiSize) {
        std::string errorMsg = "Number of MPI ranks " + std::to_string(metadata.mpiSize) + " <> "
            + std::to_string(nBoxes) + "\n";
        throw std::runtime_error(errorMsg);
    }
    netCDF::NcGroup bboxGroup(ncFile.getGroup(bboxName));
    std::vector<size_t> index(1, metadata.mpiMyRank);
    bboxGroup.getVar("global_x").getVar(index, &metadata.localCornerX);
    bboxGroup.getVar("global_y").getVar(index, &metadata.localCornerY);
    bboxGroup.getVar("local_extent_x").getVar(index, &metadata.localExtentX);
    bboxGroup.getVar("local_extent_y").getVar(index, &metadata.localExtentY);
    ncFile.close();
}
#endif

void RectGridIO::dumpModelState(const ModelState& state, const ModelMetadata& metadata,
    const std::string& filePath, bool isRestart) const
{
#ifdef USE_MPI
    auto filePathRank = filePath + "_" + std::to_string(metadata.mpiMyRank);
    netCDF::NcFile ncFile(filePathRank, netCDF::NcFile::replace);
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

    std::vector<std::string> dimensionNames = { "x", "y", "z", "t", "component", "u", "v", "w" };

    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
    netCDF::NcDim xDim = dataGroup.addDim(dimensionNames[0], nx);
    netCDF::NcDim yDim = dataGroup.addDim(dimensionNames[1], ny);
    std::vector<netCDF::NcDim> dims2 = { yDim, xDim };
    netCDF::NcDim zDim = dataGroup.addDim(dimensionNames[2], nz);
    std::vector<netCDF::NcDim> dims3 = { zDim, yDim, xDim };

    for (const auto entry : state.data) {
        const std::string& name = entry.first;
        if (entry.second.getType() == ModelArray::Type::H && entry.second.trueSize() > 0) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims2));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value);
            var.putVar(entry.second.getData());
        } else if (entry.second.getType() == ModelArray::Type::Z && entry.second.trueSize() > 0) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims3));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value);
            var.putVar(entry.second.getData());
        }
    }

    ncFile.close();
}

} /* namespace Nextsim */
