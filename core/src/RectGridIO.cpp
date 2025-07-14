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
void dimensionSetter(const netCDF::NcFile& ncfile, const std::string& fieldName,
    ModelArray::Type type, ModelMetadata& metadata)
{
    size_t nDims = ncfile.getVar(fieldName).getDimCount();
    ModelArray::MultiDim dims;
    dims.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        dims[d] = ncfile.getVar(fieldName).getDim(d).getSize();
    }
    // The dimensions in the netCDF are in the reverse order compared to ModelArray
    std::reverse(dims.begin(), dims.end());
    // Replace X, Y dimensions with local extends
    dims[0] = metadata.localExtentX;
    dims[1] = metadata.localExtentY;
    ModelArray::setDimensions(type, dims);
}
#else
void dimensionSetter(
    const netCDF::NcFile& ncfile, const std::string& fieldName, ModelArray::Type type)
{
    size_t nDims = ncfile.getVar(fieldName).getDimCount();
    ModelArray::MultiDim dims;
    dims.resize(nDims);
    for (size_t d = 0; d < nDims; ++d) {
        dims[d] = ncfile.getVar(fieldName).getDim(d).getSize();
    }
    // The dimensions in the netCDF are in the reverse order compared to ModelArray
    std::reverse(dims.begin(), dims.end());
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
    netCDF::NcFilePar ncfile(filePath, netCDF::NcFile::read, metadata.mpiComm);
#else
    netCDF::NcFile ncfile(filePath, netCDF::NcFile::read);
#endif

#ifdef USE_MPI
    // Get the sizes of the four types of field
    // HField from hice
    dimensionSetter(ncfile, hiceName, ModelArray::Type::H, metadata);
    // UField from hice
    dimensionSetter(ncfile, hiceName, ModelArray::Type::U, metadata);
    // VField from hice
    dimensionSetter(ncfile, hiceName, ModelArray::Type::V, metadata);
#else
    // Get the sizes of the four types of field
    // HField from hice
    dimensionSetter(ncfile, hiceName, ModelArray::Type::H);
    // UField from hice
    dimensionSetter(ncfile, hiceName, ModelArray::Type::U);
    // VField from hice
    dimensionSetter(ncfile, hiceName, ModelArray::Type::V);
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
    ncfile.getVar(maskName).getVar(start, size, &state.data[maskName][0]);
    state.data[hiceName] = ModelArray::HField();
    ncfile.getVar(hiceName).getVar(start, size, &state.data[hiceName][0]);
    state.data[ciceName] = ModelArray::HField();
    ncfile.getVar(ciceName).getVar(start, size, &state.data[ciceName][0]);
    state.data[hsnowName] = ModelArray::HField();
    ncfile.getVar(hsnowName).getVar(start, size, &state.data[hsnowName][0]);
    state.data[tsurfName] = ModelArray::HField();
    ncfile.getVar(tsurfName).getVar(start, size, &state.data[tsurfName][0]);
    // coordinates on the H grid
    if (ncfile.getVars().count(xName) > 0) {
        state.data[xName] = ModelArray::HField();
        ncfile.getVar(xName).getVar(start, size, &state.data[xName][0]);
        state.data[yName] = ModelArray::HField();
        ncfile.getVar(yName).getVar(start, size, &state.data[yName][0]);
    } else {
        state.data[longitudeName] = ModelArray::HField();
        ncfile.getVar(longitudeName).getVar(start, size, &state.data[longitudeName][0]);
        state.data[latitudeName] = ModelArray::HField();
        ncfile.getVar(latitudeName).getVar(start, size, &state.data[latitudeName][0]);
    }

    ncfile.close();
    return state;
}

void RectGridIO::dumpModelState(const ModelState& state, const ModelMetadata& metadata,
    const std::string& filePath, bool isRestart) const
{
#ifdef USE_MPI
    netCDF::NcFilePar ncfile(filePath, netCDF::NcFile::replace, metadata.mpiComm);
#else
    netCDF::NcFile ncfile(filePath, netCDF::NcFile::replace);
#endif

    CommonRestartMetadata::writeStructureType(ncfile, metadata);

    CommonRestartMetadata::writeRestartMetadata(ncfile, metadata);
    typedef ModelArray::Type Type;

    int nx = ModelArray::dimensions(Type::H)[0];
    int ny = ModelArray::dimensions(Type::H)[1];

    std::vector<std::string> dimensionNames = { "xdim", "ydim", "t", "component", "u", "v", "w" };

    // Create the dimension data
#ifdef USE_MPI
    netCDF::NcDim xDim = ncfile.addDim(dimensionNames[0], metadata.globalExtentX);
    netCDF::NcDim yDim = ncfile.addDim(dimensionNames[1], metadata.globalExtentY);
#else
    netCDF::NcDim xDim = ncfile.addDim(dimensionNames[0], nx);
    netCDF::NcDim yDim = ncfile.addDim(dimensionNames[1], ny);
#endif
    std::vector<netCDF::NcDim> dims2 = { yDim, xDim };
#ifdef USE_MPI
    // Set the origins and extensions for reading 2D data based
    // on MPI decomposition
    std::vector<size_t> start2 = { static_cast<size_t>(metadata.localCornerY),
        static_cast<size_t>(metadata.localCornerX) };
    std::vector<size_t> size2 = { static_cast<size_t>(metadata.localExtentY),
        static_cast<size_t>(metadata.localExtentX) };
#endif

    for (const auto entry : state.data) {
        const std::string& name = entry.first;
        if (entry.second.trueSize() > 0) {
            netCDF::NcVar var(ncfile.addVar(name, netCDF::ncDouble, dims2));
            var.putAtt(mdiName, netCDF::ncDouble, MissingData::value());
#ifdef USE_MPI
            var.putVar(start2, size2, entry.second.getData());
#else
            var.putVar(entry.second.getData());
#endif
        }
    }

    ncfile.close();
}

} /* namespace Nextsim */
