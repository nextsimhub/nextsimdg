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
static std::string sstName = "sst";
static std::string sssName = "sss";
static std::string maskName = "mask";

static const std::string mdiName = "missing_value";

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

    state[maskName] = ModelArray::HField();
    dataGroup.getVar(maskName).getVar(&state[maskName][0]);
    state[hiceName] = ModelArray::HField();
    dataGroup.getVar(hiceName).getVar(&state[hiceName][0]);
    state[ciceName] = ModelArray::HField();
    dataGroup.getVar(ciceName).getVar(&state[ciceName][0]);
    state[hsnowName] = ModelArray::HField();
    dataGroup.getVar(hsnowName).getVar(&state[hsnowName][0]);
    state[sstName] = ModelArray::HField();
    dataGroup.getVar(sstName).getVar(&state[sstName][0]);
    state[sssName] = ModelArray::HField();
    dataGroup.getVar(sssName).getVar(&state[sssName][0]);
    state[ticeName] = ModelArray::ZField();
    dataGroup.getVar(ticeName).getVar(&state[ticeName][0]);

    ncFile.close();
    return state;
}

void RectGridIO::dumpModelState(
    const ModelState& state, const ModelMetadata& metadata, const std::string& filePath, bool isRestart) const
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

    for (const auto entry : state) {
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
