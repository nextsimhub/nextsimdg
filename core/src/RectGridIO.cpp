/*!
 * @file RectGridIO.cpp
 *
 * @date Feb 8, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/RectGridIO.hpp"

#include "include/RectangularGrid.hpp"
#include "include/ElementData.hpp"
#include "include/IStructure.hpp"

#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncVar.h>

#include <map>
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

std::string hiceName = "hice";
std::string ciceName = "cice";
std::string hsnowName = "hsnow";
std::string ticeName = "tice";
std::string sstName = "sst";
std::string sssName = "sss";

typedef std::map<StringName, std::string> NameMap;

// See https://isocpp.org/wiki/faq/pointers-to-members#macro-for-ptr-to-memfn
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

// pointer to a member of PrognosticData that takes no arguments and returns a
// double. See https://isocpp.org/wiki/faq/pointers-to-members#typedef-for-ptr-to-memfn
typedef double (PrognosticData::*ProgDoubleFn)() const;

// Map between variable names and retrieval functions
// See https://isocpp.org/wiki/faq/pointers-to-members#array-memfnptrs
// clang-format off
static const std::map<std::string, ProgDoubleFn> variableFunctions
= {    { hiceName, &PrognosticData::iceThickness },
       { ciceName, &PrognosticData::iceConcentration },
       { hsnowName, &PrognosticData::snowThickness },
       { sstName, &PrognosticData::seaSurfaceTemperature },
       { sssName, &PrognosticData::seaSurfaceSalinity } };
// clang-format on

void initMeta(std::vector<ElementData>& data, RectGridIO::GridDimensions& dims, const netCDF::NcGroup& metaGroup)
{
    // No metadata to initialize
}

void initData(std::vector<ElementData>& data, RectGridIO::GridDimensions& dims, const netCDF::NcGroup& dataGroup)
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

void RectGridIO::init(std::vector<ElementData>& data, const std::string& filePath, GridDimensions& dims)
{
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);

    netCDF::NcGroup metaGroup(ncFile.getGroup(IStructure::metadataNodeName()));
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

    initMeta(data, dims, metaGroup);
    initData(data, dims, dataGroup);

    ncFile.close();
}

void dumpMeta(
    const std::vector<ElementData>& data, netCDF::NcGroup& metaGroup, const NameMap& nameMap)
{
    metaGroup.putAtt(IStructure::typeNodeName(), nameMap.at(StringName::STRUCTURE));
}

/*
 * Uses a pointer-to-member-function to access the data in PrognosticData
 * element by element. Please see the references in the comments associated
 * with ProgDoubleFn and CALL_MEMBER_FN for more details.
 */
std::vector<double> gather(const std::vector<ElementData>& data, ProgDoubleFn pFunc)
{
    std::vector<double> gathered(data.size());
    for (int i = 0; i < data.size(); ++i) {
        gathered[i] = CALL_MEMBER_FN(data[i], pFunc)();
    }
    return gathered;
}

void dumpData(const std::vector<ElementData>& data, const RectGridIO::GridDimensions& dims,
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


void RectGridIO::dump(const std::vector<ElementData>& data, const std::string& filePath, const GridDimensions& dims) const
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

} /* namespace Nextsim */
