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
    std::vector<size_t> start(2);
    std::vector<size_t> size(2);
    auto& metadata = ModelMetadata::getInstance();
    start[0] = metadata.getLocalCornerY();
    start[1] = metadata.getLocalCornerX();
    size[0] = metadata.getLocalExtentY();
    size[1] = metadata.getLocalExtentX();
#else
    std::vector<size_t> start = { 0, 0 };
    std::vector<size_t> size = ModelArray::dimensions(ModelArray::Type::H);
    std::reverse(size.begin(), size.end());
#endif
    state.data[maskName] = ModelArray::HField();
    ncFile.getVar(maskName).getVar(start, size, &state.data[maskName][0]);
    state.data[hiceName] = ModelArray::HField();
    ncFile.getVar(hiceName).getVar(start, size, &state.data[hiceName][0]);
    state.data[ciceName] = ModelArray::HField();
    ncFile.getVar(ciceName).getVar(start, size, &state.data[ciceName][0]);
    state.data[hsnowName] = ModelArray::HField();
    ncFile.getVar(hsnowName).getVar(start, size, &state.data[hsnowName][0]);
    state.data[tsurfName] = ModelArray::HField();
    ncFile.getVar(tsurfName).getVar(start, size, &state.data[tsurfName][0]);
    // coordinates on the H grid
    if (ncFile.getVars().count(xName) > 0) {
        state.data[xName] = ModelArray::HField();
        ncFile.getVar(xName).getVar(start, size, &state.data[xName][0]);
        state.data[yName] = ModelArray::HField();
        ncFile.getVar(yName).getVar(start, size, &state.data[yName][0]);
    } else {
        state.data[longitudeName] = ModelArray::HField();
        ncFile.getVar(longitudeName).getVar(start, size, &state.data[longitudeName][0]);
        state.data[latitudeName] = ModelArray::HField();
        ncFile.getVar(latitudeName).getVar(start, size, &state.data[latitudeName][0]);
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
