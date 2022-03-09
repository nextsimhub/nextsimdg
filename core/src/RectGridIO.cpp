/*!
 * @file RectGridIO.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/RectGridIO.hpp"

#include "include/ElementData.hpp"
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

static void initMeta(std::vector<ElementData>& data, RectGridIO::GridDimensions& dims,
    const netCDF::NcGroup& metaGroup)
{
    // No metadata to initialize
}

static void initData(std::vector<ElementData>& data, RectGridIO::GridDimensions& dims,
    const netCDF::NcGroup& dataGroup)
{
    // Get the number of array sizes from the dimension data of the ice temperature array
    int nDims = 3;
    int dimArray[nDims];

    for (int iDim = 0; iDim < nDims; ++iDim) {
        dimArray[iDim] = dataGroup.getVar(ticeName).getDim(iDim).getSize();
    }

    dims.nx = dimArray[0];
    dims.ny = dimArray[1];
    dims.nz = dimArray[2];

    data.resize(dims.nx * dims.ny);

    for (int i = 0; i < dims.nx; ++i) {
        for (int j = 0; j < dims.ny; ++j) {
            int linearIndex = i * dims.ny + j;
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

            std::vector<double> tice(dims.nz);
            for (int l = 0; l < dims.nz; ++l) {
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

void RectGridIO::init(
    std::vector<ElementData>& data, const std::string& filePath, GridDimensions& dims)
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);

    netCDF::NcGroup metaGroup(ncFile.getGroup(IStructure::metadataNodeName()));
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

    initMeta(data, dims, metaGroup);
    initData(data, dims, dataGroup);

    ncFile.close();
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

static void dumpData(const std::vector<ElementData>& data, const RectGridIO::GridDimensions& dims,
    netCDF::NcGroup& dataGroup, const NameMap& nameMap)
{
    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
    netCDF::NcDim xDim = dataGroup.addDim(nameMap.at(StringName::X_DIM), dims.nx);
    netCDF::NcDim yDim = dataGroup.addDim(nameMap.at(StringName::Y_DIM), dims.ny);

    std::vector<netCDF::NcDim> dims2 = { xDim, yDim };
    for (auto fnNamePair : variableFunctions) {
        const std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims2));
        std::vector<double> gathered = gather(data, variableFunctions.at(name));

        var.putVar(gathered.data());
    }

    int nLayers = dims.nz;
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

void RectGridIO::dump(const std::vector<ElementData>& data, const std::string& filePath,
    const GridDimensions& dims) const
{
    NameMap nameMap = {
        { StringName::METADATA_NODE, IStructure::metadataNodeName() },
        { StringName::DATA_NODE, IStructure::dataNodeName() },
        { StringName::STRUCTURE, RectangularGrid::structureName },
        { StringName::X_DIM, RectangularGrid::xDimName },
        { StringName::Y_DIM, RectangularGrid::yDimName },
        { StringName::Z_DIM, RectangularGrid::nIceLayersName },
    };
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);

    netCDF::NcGroup metaGroup = ncFile.addGroup(nameMap.at(StringName::METADATA_NODE));
    netCDF::NcGroup dataGroup = ncFile.addGroup(nameMap.at(StringName::DATA_NODE));

    dumpMeta(data, metaGroup, nameMap);
    dumpData(data, dims, dataGroup, nameMap);

    ncFile.close();
}

void dimensionSetter(const netCDF::NcGroup& dataGroup, const std::string& fieldName, ModelArray::Type type)
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
            { Type::H, {hiceName, ciceName, hsnowName, sstName, sssName} },
            { Type::U, { } },
            { Type::V, { } },
            { Type::Z, {ticeName} },
    };

    for (auto entry : fields) {
        Type type = entry.first;
        // Create the dimension data
        std::vector<netCDF::NcDim> dimVec;
        for (size_t d = 0; d < ModelArray::nDimensions(type); ++d) {

        }
    }
//    dumpData(data, dims, dataGroup, nameMap);

    ncFile.close();

}


} /* namespace Nextsim */
