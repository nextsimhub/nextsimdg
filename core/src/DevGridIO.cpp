/*!
 * @file DevGridIO.cpp
 *
 * @date Jan 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevGridIO.hpp"

#include "include/DevGrid.hpp"
#include "include/IStructure.hpp"
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
static std::string sstName = "sst";
static std::string sssName = "sss";

typedef std::map<StringName, std::string> NameMap;

static const std::string metaName = "meta";
static const std::string dataName = "data";

// Metadata initialization
static void initModelMetaData(const netCDF::NcGroup& metaGroup) { }

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
            HField data = ModelArray::HField(varName);
            data.setData(buffer.data());
            auto [i, y] = state.insert({ varName, data });
        } else if (nDims == 3) {
            ZField data = ModelArray::ZField(varName);
            data.setData(buffer.data());
            state[varName] = data;
        }
    }

    return state;
}

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

void dumpModelMeta(const ModelState& state, netCDF::NcGroup& metaGroup)
{
    metaGroup.putAtt(IStructure::typeNodeName(), DevGrid::structureName);
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

    for (const auto entry : state) {
        const std::string& name = entry.first;
        if (entry.second.getType() == ModelArray::Type::H) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims2));
            var.putVar(entry.second.getData());
        } else if (entry.second.getType() == ModelArray::Type::Z) {
            netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims3));
            var.putVar(entry.second.getData());
        }
    }
}

void DevGridIO::dumpModelState(const ModelState& state, const std::string& filePath) const
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);
    netCDF::NcGroup metaGroup = ncFile.addGroup(metaName);
    netCDF::NcGroup dataGroup = ncFile.addGroup(dataName);
    dumpModelMeta(state, metaGroup);
    dumpModelData(state, dataGroup);
    ncFile.close();
}

} /* namespace Nextsim */
