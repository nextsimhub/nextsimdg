/*!
 * @file DevGridIO.cpp
 *
 * @date Jan 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Athena Elafrou <ae488@cam.ac.uk>
 */

#include "include/DevGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/DevGrid.hpp"
#include "include/IStructure.hpp"
#include "include/MissingData.hpp"
#include "include/ModelArray.hpp"

#include <cstddef>
#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncVar.h>

#include <map>
#include <string>
#include <vector>

namespace Nextsim {

// Forward declarations
enum class StringName {
    METADATA_NODE,
    DATA_NODE,
    STRUCTURE,
    X_DIM,
    Y_DIM,
    Z_DIM,
};

static std::string hiceName = "hice";
static std::string ciceName = "cice";
static std::string hsnowName = "hsnow";
static std::string ticeName = "tice";

typedef std::map<StringName, std::string> NameMap;

static const std::string metaName = "meta";
static const std::string dataName = "data";
static const std::string mdiName = "missing_value";
#ifdef USE_MPI
static const std::string bboxName = "bounding_boxes";
#endif // USE_MPI

// Metadata initialization
static void initModelMetaData(const netCDF::NcGroup& metaGroup) { }

#ifdef USE_MPI
static ModelState initModelData(
    const netCDF::NcGroup& dataGroup, const netCDF::NcGroup& bboxGroup, MPI_Comm mpiComm)
{
    // Get the global array sizes from the dimension data of the ice temperature array
    size_t nDims = 3;
    std::vector<size_t> globalDim(nDims);
    for (int iDim = 0; iDim < nDims; ++iDim) {
        globalDim[iDim] = dataGroup.getVar(ticeName).getDim(iDim).getSize();
    }

    // Get the local array sizes for this process
    int nProcs = bboxGroup.getDim("P").getSize();
    int mpiSize, mpiRank;
    MPI_Comm_size(mpiComm, &mpiSize);
    MPI_Comm_rank(mpiComm, &mpiRank);
    assert(mpiSize == nProcs);
    std::vector<size_t> index(1, mpiRank);
    int topX, topY, cntX, cntY;
    bboxGroup.getVar("global_x").getVar(index, &topX);
    bboxGroup.getVar("global_y").getVar(index, &topY);
    bboxGroup.getVar("local_extent_x").getVar(index, &cntX);
    bboxGroup.getVar("local_extent_y").getVar(index, &cntY);

    std::vector<size_t> localDim3(nDims);
    std::vector<size_t> localDim2(nDims - 1);
    localDim3[0] = localDim2[0] = cntX;
    localDim3[1] = localDim2[1] = cntY;
    localDim3[2] = globalDim[2];
    // assert dimensions match

    ModelArray::setDimensions(ModelArray::Type::H, localDim2);
    ModelArray::setDimensions(ModelArray::Type::U, localDim2);
    ModelArray::setDimensions(ModelArray::Type::V, localDim2);
    ModelArray::setDimensions(ModelArray::Type::Z, localDim3);

    // Setup information to load my part of the model data
    std::vector<size_t> start(2);
    std::vector<size_t> count(2);
    // Coordinate of first element
    start[0] = topX;
    start[1] = topY;
    // Number of elements in every dimension
    count[0] = cntX;
    count[1] = cntY;

    // Loop over all data fields, add their name and contents to the ModelState
    std::multimap<std::string, netCDF::NcVar> varMap = dataGroup.getVars();
    ModelState state;
    for (const auto var : varMap) {
        std::string varName = var.first;
        std::vector<double> buffer(cntX * cntY);
        var.second.getVar(start, count, buffer.data());
        int nDims = var.second.getDimCount();
        if (nDims == 2) {
            HField data = ModelArray::HField();
            data.setData(buffer.data());
            auto [i, y] = state.data.insert({ varName, data });
        } else if (nDims == 3) {
            ZField data = ModelArray::ZField();
            data.setData(buffer.data());
            state.data[varName] = data;
        }
    }

    return state;
}
#else
static ModelState initModelData(const netCDF::NcGroup& dataGroup)
{
    // Get the number of array sizes from the dimension data of the ice temperature array
    size_t nDims = 3;
    std::vector<size_t> dim3(nDims);
    std::vector<size_t> dim2(nDims - 1);

    for (int iDim = 0; iDim < nDims; ++iDim) {
        dim3[iDim] = dataGroup.getVar(ticeName).getDim(iDim).getSize();
    }
    dim2[0] = dim3[0];
    dim2[1] = dim3[1];

    ModelArray::setDimensions(ModelArray::Type::H, dim2);
    ModelArray::setDimensions(ModelArray::Type::U, dim2);
    ModelArray::setDimensions(ModelArray::Type::V, dim2);
    ModelArray::setDimensions(ModelArray::Type::Z, dim3);

    // Loop over all data fields, add their name and contents to the ModelState
    std::multimap<std::string, netCDF::NcVar> varMap = dataGroup.getVars();
    ModelState state;
    for (const auto var : varMap) {
        // Get the dimensions and read into an intermediate buffer
        int nDims = var.second.getDimCount();
        std::vector<netCDF::NcDim> dims = var.second.getDims();
        size_t totalSz = 1;
        for (const auto dim : dims) {
            totalSz *= dim.getSize();
        }

        std::string varName = var.first;
        std::vector<double> buffer(totalSz);
        var.second.getVar(buffer.data());
        if (nDims == 2) {
            HField data = ModelArray::HField();
            data.setData(buffer.data());
            auto [i, y] = state.data.insert({ varName, data });
        } else if (nDims == 3) {
            ZField data = ModelArray::ZField();
            data.setData(buffer.data());
            state.data[varName] = data;
        }
    }

    return state;
}
#endif // USE_MPI

#ifdef USE_MPI
ModelState DevGridIO::getModelState(
    const std::string& modelFilePath, const std::string& partitionFilePath) const
{
    netCDF::NcFile modelNcFile(modelFilePath, netCDF::NcFile::read);
    netCDF::NcGroup metaGroup(modelNcFile.getGroup(metaName));
    initModelMetaData(metaGroup);

    netCDF::NcGroup dataGroup(modelNcFile.getGroup(dataName));
    netCDF::NcFile partitionNcFile(partitionFilePath, netCDF::NcFile::read);
    netCDF::NcGroup bboxGroup(partitionNcFile.getGroup(bboxName));
    ModelState ms = initModelData(dataGroup, bboxGroup, grid->getComm());

    modelNcFile.close();
    partitionNcFile.close();

    return ms;
}
#else
ModelState DevGridIO::getModelState(const std::string& filePath) const
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);

    netCDF::NcGroup metaGroup(ncFile.getGroup(metaName));
    netCDF::NcGroup dataGroup(ncFile.getGroup(dataName));
    initModelMetaData(metaGroup);
    ModelState ms = initModelData(dataGroup);

    ncFile.close();

    return ms;
}
#endif // USE_MPI

void dumpModelMeta(const ModelMetadata& metadata, netCDF::NcGroup& metaGroup)
{
    CommonRestartMetadata::writeRestartMetadata(metaGroup, metadata);
}

void dumpModelData(const ModelState& state, netCDF::NcGroup& dataGroup)
{
    int nx = DevGrid::nx;
    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
    netCDF::NcDim xDim = dataGroup.addDim(DevGrid::xDimName, nx);
    netCDF::NcDim yDim = dataGroup.addDim(DevGrid::yDimName, nx);
    std::vector<netCDF::NcDim> dims2 = { xDim, yDim };
    int nLayers = 1;
    netCDF::NcDim zDim = dataGroup.addDim(DevGrid::nIceLayersName, nLayers);
    std::vector<netCDF::NcDim> dims3 = { xDim, yDim, zDim };

    for (const auto entry : state.data) {
        const std::string& name = entry.first;
        if (entry.second.getType() == ModelArray::Type::H) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims2));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
            var.putVar(entry.second.getData());
        } else if (entry.second.getType() == ModelArray::Type::Z) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims3));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
            var.putVar(entry.second.getData());
        }
    }
}

void DevGridIO::dumpModelState(const ModelState& state, const ModelMetadata& metadata,
    const std::string& filePath, bool isRestart) const
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);
    netCDF::NcGroup metaGroup = ncFile.addGroup(metaName);
    netCDF::NcGroup dataGroup = ncFile.addGroup(dataName);
    CommonRestartMetadata::writeStructureType(ncFile, metadata);
    dumpModelMeta(metadata, metaGroup);
    dumpModelData(state, dataGroup);
    ncFile.close();
}

} /* namespace Nextsim */
