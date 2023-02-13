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

#ifdef USE_MPI
ModelState RectGridIO::getModelState(const std::string& filePath, const std::string& partitionFile)
#else
ModelState RectGridIO::getModelState(const std::string& filePath)
#endif
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
    // Since the ZFierld might not have the same dimensions as the tice field
    // in the file, a little more work is required.
    state.data[ticeName] = ModelArray::ZField();
    std::vector<size_t> startVector = { 0, 0, 0 };
    std::vector<size_t> zArrayDims = ModelArray::dimensions(ModelArray::Type::Z);
    std::reverse(zArrayDims.begin(), zArrayDims.end());
    dataGroup.getVar(ticeName).getVar(startVector, zArrayDims, &state.data[ticeName][0]);

    ncFile.close();
    return state;
}

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
