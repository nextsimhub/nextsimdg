/*!
 * @file RectGridIO.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/RectGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/IStructure.hpp"
#include "include/MissingData.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelState.hpp"
#include "include/RectangularGrid.hpp"

#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncVar.h>

#include <list>
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
static std::string maskName = "mask";

static const std::string mdiName = "missing_value";
#ifdef USE_MPI
static const std::string bboxName = "bounding_boxes";
#endif // USE_MPI

typedef std::map<StringName, std::string> NameMap;

void dimensionSetter(
    const netCDF::NcGroup& dataGroup, const std::string& fieldName, ModelArray::Type type)
{
    size_t nDims = dataGroup.getVar(fieldName).getDimCount();
    ModelArray::Dimensions dims;
    dims.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        dims[d] = dataGroup.getVar(fieldName).getDim(d).getSize();
    }
    ModelArray::setDimensions(type, dims);
}

#ifdef USE_MPI
ModelState RectGridIO::getModelState(
    const std::string& modelFilePath, const std::string& partitionFilePath) const
{
    ModelState state;
    netCDF::NcFile modelNcFile(modelFilePath, netCDF::NcFile::read);
    netCDF::NcGroup modelGroup(modelNcFile.getGroup(IStructure::dataNodeName()));

    // Get the global array sizes from the dimension data of the ice temperature array
    // We are assuming that H,U,V type arrays share the same dimensions and Z type arrays share the
    // first 2 dimesions, but have an extra dimension for the number of layers
    size_t nDims = modelGroup.getVar(ticeName).getDimCount();
    ModelArray::Dimensions globalDim;
    globalDim.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        globalDim[d] = modelGroup.getVar(ticeName).getDim(d).getSize();
    }

    // Get the local sizes of the four types of field
    netCDF::NcFile partitionNcFile(partitionFilePath, netCDF::NcFile::read);
    netCDF::NcGroup bboxGroup(partitionNcFile.getGroup(bboxName));
    int nProcs = bboxGroup.getDim("P").getSize();
    MPI_Comm mpiComm = grid->getComm();
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
    partitionNcFile.close();

    std::vector<size_t> localDim3(nDims);
    std::vector<size_t> localDim2(nDims - 1);
    localDim3[0] = localDim2[0] = cntX;
    localDim3[1] = localDim2[1] = cntY;
    localDim3[2] = globalDim[2];
    // assert dimensions match

    // HField from hice
    ModelArray::setDimensions(ModelArray::Type::H, localDim2);
    // UField from hice        TODO replace with u velocity once it is present
    ModelArray::setDimensions(ModelArray::Type::U, localDim2);
    // VField from hice        TODO replace with v velocity once it is present
    ModelArray::setDimensions(ModelArray::Type::V, localDim2);
    // ZField from tice
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

    state.data[maskName] = ModelArray::HField();
    modelGroup.getVar(maskName).getVar(start, count, &state.data[maskName][0]);
    state.data[hiceName] = ModelArray::HField();
    modelGroup.getVar(hiceName).getVar(start, count, &state.data[hiceName][0]);
    state.data[ciceName] = ModelArray::HField();
    modelGroup.getVar(ciceName).getVar(start, count, &state.data[ciceName][0]);
    state.data[hsnowName] = ModelArray::HField();
    modelGroup.getVar(hsnowName).getVar(start, count, &state.data[hsnowName][0]);
    state.data[ticeName] = ModelArray::ZField();
    // Set third dimension for Z type array
    start.push_back(0);
    count.push_back(localDim3[2]);
    modelGroup.getVar(ticeName).getVar(start, count, &state.data[ticeName][0]);

    modelNcFile.close();

    // FIXME: Also update RectagularGrid's private data members

    return state;
}
#else
ModelState RectGridIO::getModelState(const std::string& filePath)
{
    ModelState state;
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

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
    state.data[ticeName] = ModelArray::ZField();
    dataGroup.getVar(ticeName).getVar(&state.data[ticeName][0]);

    ncFile.close();

    return state;
}
#endif // USE_MPI

void RectGridIO::dumpModelState(const ModelState& state, const ModelMetadata& metadata,
    const std::string& filePath, bool isRestart) const
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);

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
    std::vector<netCDF::NcDim> dims2 = { xDim, yDim };
    netCDF::NcDim zDim = dataGroup.addDim(dimensionNames[2], nz);
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

    ncFile.close();
}

} /* namespace Nextsim */
