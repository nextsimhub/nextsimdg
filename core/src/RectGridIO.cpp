/*!
 * @file RectGridIO.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/RectGridIO.hpp"

#include "include/IStructure.hpp"
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

    state[hiceName] = ModelArray::HField(hiceName);
    dataGroup.getVar(hiceName).getVar(&state[hiceName][0]);
    state[ciceName] = ModelArray::HField(ciceName);
    dataGroup.getVar(ciceName).getVar(&state[ciceName][0]);
    state[hsnowName] = ModelArray::HField(hsnowName);
    dataGroup.getVar(hsnowName).getVar(&state[hsnowName][0]);
    state[sstName] = ModelArray::HField(sstName);
    dataGroup.getVar(sstName).getVar(&state[sstName][0]);
    state[sssName] = ModelArray::HField(sssName);
    dataGroup.getVar(sssName).getVar(&state[sssName][0]);
    state[ticeName] = ModelArray::ZField(ticeName);
    dataGroup.getVar(ticeName).getVar(&state[ticeName][0]);

    ncFile.close();
    return state;
}

void RectGridIO::dumpModelState(const ModelState& state, const std::string& filePath) const
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);

    netCDF::NcGroup metaGroup = ncFile.addGroup(IStructure::metadataNodeName());
    netCDF::NcGroup dataGroup = ncFile.addGroup(IStructure::dataNodeName());

    metaGroup.putAtt(IStructure::typeNodeName(), RectangularGrid::structureName);

    typedef ModelArray::Type Type;

    std::map<Type, std::list<std::string>> fields = {
        { Type::H, { hiceName, ciceName, hsnowName, sstName, sssName } },
        { Type::U, {} },
        { Type::V, {} },
        { Type::Z, { ticeName } },
    };

    std::vector<std::string> dimensionNames = { "x", "y", "z", "t", "component", "u", "v", "w" };

    for (auto entry : fields) {
        if (entry.second.size() == 0)
            continue;
        Type type = entry.first;

        // Create the dimension data
        std::vector<netCDF::NcDim> dimVec;
        for (size_t d = 0; d < ModelArray::nDimensions(type); ++d) {
            dimVec.push_back(dataGroup.addDim(ModelArray::typeNames.at(type) + dimensionNames[d],
                ModelArray::dimensions(type)[d]));
        }

        // Write out the data for each field of this type
        for (auto field : entry.second) {
            netCDF::NcVar var(dataGroup.addVar(field, netCDF::ncDouble, dimVec));
            var.putVar(&state.at(field)[0]);
        }
    }
    //    dumpData(data, dims, dataGroup, nameMap);

    ncFile.close();
}

} /* namespace Nextsim */
