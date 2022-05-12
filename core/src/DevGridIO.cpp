/*!
 * @file DevGridIO.cpp
 *
 * @date Jan 14, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevGridIO.hpp"

#include "include/DevGrid.hpp"
#include "include/ElementData.hpp"
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

static void initGroup(std::vector<ElementData>& data, netCDF::NcGroup& grp, const NameMap& nameMap);
static void dumpGroup(
    const std::vector<ElementData>& data, netCDF::NcGroup& grp, const NameMap& nameMap);

// See https://isocpp.org/wiki/faq/pointers-to-members#macro-for-ptr-to-memfn
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

// pointer to a member of PrognosticElementData that takes no arguments and returns a
// double. See https://isocpp.org/wiki/faq/pointers-to-members#typedef-for-ptr-to-memfn
typedef double (PrognosticElementData::*ProgDoubleFn)() const;

// Map between variable names and retrieval functions
// See https://isocpp.org/wiki/faq/pointers-to-members#array-memfnptrs
// clang-format off
static const std::map<std::string, ProgDoubleFn> variableFunctions
= {    { hiceName, &PrognosticElementData::iceThickness },
       { ciceName, &PrognosticElementData::iceConcentration },
       { hsnowName, &PrognosticElementData::snowThickness },
       { sstName, &PrognosticElementData::seaSurfaceTemperature },
       { sssName, &PrognosticElementData::seaSurfaceSalinity } };
// clang-format on

void DevGridIO::init(std::vector<ElementData>& data, const std::string& filePath) const
{
    NameMap nameMap = {
        { StringName::METADATA_NODE, IStructure::metadataNodeName() },
        { StringName::DATA_NODE, IStructure::dataNodeName() },
        { StringName::STRUCTURE, DevGrid::structureName },
        { StringName::X_DIM, DevGrid::xDimName },
        { StringName::Y_DIM, DevGrid::yDimName },
        { StringName::Z_DIM, DevGrid::nIceLayersName },
    };
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
    initGroup(data, ncFile, nameMap);
    ncFile.close();
}

void DevGridIO::dump(const std::vector<ElementData>& data, const std::string& filePath) const
{
    NameMap nameMap = {
        { StringName::METADATA_NODE, IStructure::metadataNodeName() },
        { StringName::DATA_NODE, IStructure::dataNodeName() },
        { StringName::STRUCTURE, DevGrid::structureName },
        { StringName::X_DIM, DevGrid::xDimName },
        { StringName::Y_DIM, DevGrid::yDimName },
        { StringName::Z_DIM, DevGrid::nIceLayersName },
    };
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);
    dumpGroup(data, ncFile, nameMap);
    ncFile.close();
}

static void initMeta(std::vector<ElementData>& data, const netCDF::NcGroup& metaGroup)
{
    int nx = DevGrid::nx;
    data.resize(nx * nx);
}

static void initData(std::vector<ElementData>& data, const netCDF::NcGroup& dataGroup)
{
    // Get the number of ice layers from the ice temperature data
    const int layersDim = 2;
    int nLayers = dataGroup.getVar(ticeName).getDim(layersDim).getSize();
    int nx = DevGrid::nx;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < nx; ++j) {
            int linearIndex = i * nx + j;
            double hice;
            double cice;
            double hsnow;
            double sst;
            double sss;
            std::vector<std::size_t> loc = { std::size_t(i), std::size_t(j) };
            dataGroup.getVar(hiceName).getVar(loc, &hice);
            dataGroup.getVar(ciceName).getVar(loc, &cice);
            dataGroup.getVar(hsnowName).getVar(loc, &hsnow);
            dataGroup.getVar(sstName).getVar(loc, &sst);
            dataGroup.getVar(sssName).getVar(loc, &sss);
            // Retrieve ice temperature data

            std::vector<double> tice(nLayers);
            for (int l = 0; l < nLayers; ++l) {
                std::vector<std::size_t> loc3 = { loc[0], loc[1], std::size_t(l) };
                dataGroup.getVar(ticeName).getVar(loc3, &tice[l]);
            }
            data[linearIndex]
                = PrognosticGenerator().hice(hice).cice(cice).sst(sst).sss(sss).hsnow(hsnow).tice(
                    tice);
        }
    }

    for (auto fnNamePair : variableFunctions) {
        const std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.getVar(name));
    }
}

static void initGroup(std::vector<ElementData>& data, netCDF::NcGroup& grp, const NameMap& nameMap)
{
    netCDF::NcGroup metaGroup(grp.getGroup(nameMap.at(StringName::METADATA_NODE)));
    netCDF::NcGroup dataGroup(grp.getGroup(nameMap.at(StringName::DATA_NODE)));

    initMeta(data, metaGroup);
    initData(data, dataGroup);
}

static void dumpMeta(
    const std::vector<ElementData>& data, netCDF::NcGroup& metaGroup, const NameMap& nameMap)
{
    metaGroup.putAtt(IStructure::typeNodeName(), nameMap.at(StringName::STRUCTURE));
}

/*
 * Uses a pointer-to-member-function to access the data in PrognosticElementData
 * element by element. Please see the references in the comments associated
 * with ProgDoubleFn and CALL_MEMBER_FN for more details.
 */
static std::vector<double> gather(const std::vector<ElementData>& data, ProgDoubleFn pFunc)
{
    std::vector<double> gathered(data.size());
    for (int i = 0; i < data.size(); ++i) {
        gathered[i] = CALL_MEMBER_FN(data[i], pFunc)();
    }
    return gathered;
}

static void dumpData(
    const std::vector<ElementData>& data, netCDF::NcGroup& dataGroup, const NameMap& nameMap)
{
    int nx = DevGrid::nx;
    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
    netCDF::NcDim xDim = dataGroup.addDim(nameMap.at(StringName::X_DIM), nx);
    netCDF::NcDim yDim = dataGroup.addDim(nameMap.at(StringName::Y_DIM), nx);

    std::vector<netCDF::NcDim> dims2 = { xDim, yDim };
    for (auto fnNamePair : variableFunctions) {
        const std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims2));
        std::vector<double> gathered = gather(data, variableFunctions.at(name));
        var.putVar(gathered.data());
    }

    int nLayers = data[0].nIceLayers();
    netCDF::NcDim zDim = dataGroup.addDim(nameMap.at(StringName::Z_DIM), nLayers);
    std::vector<netCDF::NcDim> dims3 = { xDim, yDim, zDim };
    netCDF::NcVar iceT(dataGroup.addVar(ticeName, netCDF::ncDouble, dims3));

    // Gather the three dimensional data explicitly (until there is more than
    // one three dimensional dataset).
    std::vector<double> tice(data.size() * nLayers);
    for (int i = 0; i < data.size(); ++i) {
        int ii = nLayers * i;
        for (int l = 0; l < nLayers; ++l) {
            tice[ii + l] = data[i].iceTemperature(l);
        }
    }
    iceT.putVar(tice.data());
}

static void dumpGroup(
    const std::vector<ElementData>& data, netCDF::NcGroup& headGroup, const NameMap& nameMap)
{
    netCDF::NcGroup metaGroup = headGroup.addGroup(nameMap.at(StringName::METADATA_NODE));
    netCDF::NcGroup dataGroup = headGroup.addGroup(nameMap.at(StringName::DATA_NODE));
    dumpMeta(data, metaGroup, nameMap);
    dumpData(data, dataGroup, nameMap);
}

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
