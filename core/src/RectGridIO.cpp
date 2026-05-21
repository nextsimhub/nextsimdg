/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Kacper Kornet <kk562@cam.ac.uk>
 */

#include "include/RectGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/IStructure.hpp"
#include "include/MissingData.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelMPI.hpp"
#include "include/ModelState.hpp"
#include "include/NetCDFUtils.hpp"
#include "include/gridNames.hpp"

#ifdef USE_MPI
#include "include/ParallelNetcdfFile.hpp"
#endif

#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncFloat.h>
#include <ncVar.h>

#include <algorithm>
#include <list>
#include <map>
#include <string>
#include <vector>

namespace Nextsim {

#ifdef USE_MPI
void dimensionSetter(
    const netCDF::NcFile& ncFile, const std::string& fieldName, ModelArray::Type type)
{
    size_t nDims = ncFile.getVar(fieldName).getDimCount();
    ModelArray::MultiDim dims;
    dims.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        dims[d] = ncFile.getVar(fieldName).getDim(d).getSize();
    }
    // The dimensions in the netCDF are in the reverse order compared to ModelArray
    std::reverse(dims.begin(), dims.end());
    // Replace X, Y dimensions with local extends
    auto& metadata = ModelMetadata::getInstance();
    dims[0] = metadata.getLocalExtentX();
    dims[1] = metadata.getLocalExtentY();
    ModelArray::setDimensions(type, dims);
}
#else
void dimensionSetter(
    const netCDF::NcFile& ncFile, const std::string& fieldName, ModelArray::Type type)
{
    size_t nDims = ncFile.getVar(fieldName).getDimCount();
    ModelArray::MultiDim dims;
    dims.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        dims[d] = ncFile.getVar(fieldName).getDim(d).getSize();
    }
    // The dimensions in the netCDF are in the reverse order compared to ModelArray
    std::reverse(dims.begin(), dims.end());
    ModelArray::setDimensions(type, dims);
}
#endif

ModelState RectGridIO::getModelState(const std::string& filePath)
{
    ModelState state;
#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance();
    netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::read, modelMPI.getComm());
#else
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
#endif

    // Get the sizes of the three types of field
    // HField from hice
    dimensionSetter(ncFile, hiceName, ModelArray::Type::H);
    // UField from hice
    dimensionSetter(ncFile, hiceName, ModelArray::Type::U);
    // VField from hice
    dimensionSetter(ncFile, hiceName, ModelArray::Type::V);

#ifdef USE_MPI
    // Set the origins and extensions for reading 2D data based
    // on MPI decomposition
    const auto& metadata = ModelMetadata::getInstance();
    const std::vector<size_t> start = { metadata.getLocalCornerY(), metadata.getLocalCornerX() };
    const std::vector<size_t> size = { metadata.getLocalExtentY(), metadata.getLocalExtentX() };
#else
    const std::vector<size_t> start = { 0, 0 };
    const std::vector<size_t> size(ModelArray::dimensions(ModelArray::Type::H).rbegin(),
        ModelArray::dimensions(ModelArray::Type::H).rend());
#endif
    state.data[maskName] = ModelArray::HField();
    readNetCDFVar(ncFile.getVar(maskName), start, size, &state.data[maskName][0]);
    state.data[hiceName] = ModelArray::HField();
    readNetCDFVar(ncFile.getVar(hiceName), start, size, &state.data[hiceName][0]);
    state.data[ciceName] = ModelArray::HField();
    readNetCDFVar(ncFile.getVar(ciceName), start, size, &state.data[ciceName][0]);
    state.data[hsnowName] = ModelArray::HField();
    readNetCDFVar(ncFile.getVar(hsnowName), start, size, &state.data[hsnowName][0]);
    state.data[tsurfName] = ModelArray::HField();
    readNetCDFVar(ncFile.getVar(tsurfName), start, size, &state.data[tsurfName][0]);
    // coordinates on the H grid
    if (ncFile.getVars().count(xName) > 0) {
        state.data[xName] = ModelArray::HField();
        readNetCDFVar(ncFile.getVar(xName), start, size, &state.data[xName][0]);
        state.data[yName] = ModelArray::HField();
        readNetCDFVar(ncFile.getVar(yName), start, size, &state.data[yName][0]);
    } else {
        state.data[longitudeName] = ModelArray::HField();
        readNetCDFVar(ncFile.getVar(longitudeName), start, size, &state.data[longitudeName][0]);
        state.data[latitudeName] = ModelArray::HField();
        readNetCDFVar(ncFile.getVar(latitudeName), start, size, &state.data[latitudeName][0]);
    }

    ncFile.close();
    return state;
}

void RectGridIO::dumpModelState(
    const ModelState& state, const std::string& filePath, bool isRestart) const
{
#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance();
    netCDF::NcFilePar ncFile(filePath, netCDF::NcFile::replace, modelMPI.getComm());
#else
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);
#endif

    CommonRestartMetadata::writeStructureType(ncFile);
    CommonRestartMetadata::writeRestartMetadata(ncFile);
    typedef ModelArray::Type Type;

    int nx = ModelArray::dimensions(Type::H)[0];
    int ny = ModelArray::dimensions(Type::H)[1];

    std::vector<std::string> dimensionNames = { "xdim", "ydim", "t", "component", "u", "v", "w" };

    // Create the dimension data
#ifdef USE_MPI
    auto& metadata = ModelMetadata::getInstance();
    netCDF::NcDim xDim = ncFile.addDim(dimensionNames[0], metadata.getGlobalExtentX());
    netCDF::NcDim yDim = ncFile.addDim(dimensionNames[1], metadata.getGlobalExtentY());
#else
    netCDF::NcDim xDim = ncFile.addDim(dimensionNames[0], nx);
    netCDF::NcDim yDim = ncFile.addDim(dimensionNames[1], ny);
#endif
    std::vector<netCDF::NcDim> dims2 = { yDim, xDim };
#ifdef USE_MPI
    // Set the origins and extensions for reading 2D data based
    // on MPI decomposition
    std::vector<size_t> start2 = { static_cast<size_t>(metadata.getLocalCornerY()),
        static_cast<size_t>(metadata.getLocalCornerX()) };
    std::vector<size_t> size2 = { static_cast<size_t>(metadata.getLocalExtentY()),
        static_cast<size_t>(metadata.getLocalExtentX()) };
#endif

    for (const auto entry : state.data) {
        const std::string& name = entry.first;
        if (entry.second.trueSize() > 0) {
            netCDF::NcVar var(ncFile.addVar(name, ToNetCDFType<FloatType>::get(), dims2));
            var.putAtt(mdiName, ToNetCDFType<FloatType>::get(), MissingData::value());
#ifdef USE_MPI
            var.putVar(start2, size2, entry.second.getData());
#else
            var.putVar(entry.second.getData());
#endif
        }
    }

    ncFile.close();
}

} /* namespace Nextsim */
