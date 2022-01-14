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

#include <cstddef>
#include <ncDim.h>
#include <ncDouble.h>
#include <ncFile.h>
#include <ncVar.h>

#include "/opt/home/include/ncDim.h"
#include "/opt/home/include/ncDouble.h"
#include "/opt/home/include/ncFile.h"
#include "/opt/home/include/ncVar.h"

#include <map>
#include <vector>

namespace Nextsim {

// Forward declarations of functions that are used herein
enum class StringName {
    METADATA_NODE,
    DATA_NODE,
    STRUCTURE,
    X_DIM,
    Y_DIM,
};

typedef std::map<StringName, std::string> NameMap;

void initGroup(std::vector<ElementData>& data, netCDF::NcGroup& grp, const NameMap& nameMap);
void dumpGroup(const std::vector<ElementData>& data, netCDF::NcGroup& grp, const NameMap& nameMap);

// See https://isocpp.org/wiki/faq/pointers-to-members#macro-for-ptr-to-memfn
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

// pointer to a member of PrognosticData that takes no arguments and returns a
// double. See https://isocpp.org/wiki/faq/pointers-to-members#typedef-for-ptr-to-memfn
typedef double (PrognosticData::*ProgDoubleFn)() const;

// Map between variable names and retrieval functions
// See https://isocpp.org/wiki/faq/pointers-to-members#array-memfnptrs
// clang-format off
static const std::map<std::string, ProgDoubleFn> variableFunctions
= {    { "hice", &PrognosticData::iceThickness },
       { "cice", &PrognosticData::iceConcentration },
       { "hsnow", &PrognosticData::snowThickness },
       { "sst", &PrognosticData::seaSurfaceTemperature },
       { "sss", &PrognosticData::seaSurfaceSalinity } };
// clang-format on


void DevGridIO::init(std::vector<ElementData>& data, const std::string& filePath) const
{
    NameMap nameMap = { {StringName::METADATA_NODE, grid->metadataNodeName},
            {StringName::DATA_NODE, grid->dataNodeName},
            {StringName::STRUCTURE, DevGrid::ourStructureName},
            {StringName::X_DIM, DevGrid::xDimName},
            {StringName::Y_DIM, DevGrid::yDimName},
    };
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::read);
    initGroup(data, ncFile, nameMap);
    ncFile.close();
}

void DevGridIO::dump(const std::vector<ElementData>& data, const std::string& filePath) const
{
    NameMap nameMap = { {StringName::METADATA_NODE, grid->metadataNodeName},
            {StringName::DATA_NODE, grid->dataNodeName},
            {StringName::STRUCTURE, DevGrid::ourStructureName},
            {StringName::X_DIM, DevGrid::xDimName},
            {StringName::Y_DIM, DevGrid::yDimName},
    };
    netCDF::NcFile ncFile(filePath, netCDF::NcFile::replace);
    dumpGroup(data, ncFile, nameMap);
    ncFile.close();
}

void initMeta(std::vector<ElementData>& data, const netCDF::NcGroup& metaGroup)
{
    int nx = DevGrid::nx;
    data.resize(nx * nx);
}

void initData(std::vector<ElementData>& data, const netCDF::NcGroup& dataGroup)
{

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
            dataGroup.getVar("hice").getVar(loc, &hice);
            dataGroup.getVar("cice").getVar(loc, &cice);
            dataGroup.getVar("hsnow").getVar(loc, &hsnow);
            dataGroup.getVar("sst").getVar(loc, &sst);
            dataGroup.getVar("sss").getVar(loc, &sss);
            // TODO How to store ice temperature data?
            std::array<double, N_ICE_TEMPERATURES> tice = { 0., 0., 0. };
            data[linearIndex] = PrognosticData::generate(hice, cice, sst, sss, hsnow, tice);
        }
    }

    for (auto fnNamePair : variableFunctions) {
        const std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.getVar(name));
    }
}

void initGroup(std::vector<ElementData>& data, netCDF::NcGroup& grp, const NameMap& nameMap)
{
    netCDF::NcGroup metaGroup(grp.getGroup(nameMap.at(StringName::METADATA_NODE)));
    netCDF::NcGroup dataGroup(grp.getGroup(nameMap.at(StringName::DATA_NODE)));

    initMeta(data, metaGroup);
    initData(data, dataGroup);
}

void dumpMeta(const std::vector<ElementData>& data, netCDF::NcGroup& metaGroup, const NameMap& nameMap)
{
    metaGroup.putAtt(IStructure::typeNodeName(), nameMap.at(StringName::STRUCTURE));
}

std::vector<double> gather(const std::vector<ElementData>& data, ProgDoubleFn pFunc)
{
    std::vector<double> gathered(data.size());
    for (int i = 0; i < data.size(); ++i) {
        gathered[i] = CALL_MEMBER_FN(data[i], pFunc)();
    }
    return gathered;
}

void dumpData(const std::vector<ElementData>& data, netCDF::NcGroup& dataGroup, const NameMap& nameMap)
{
    int nx = DevGrid::nx;
    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
    netCDF::NcDim xDim = dataGroup.addDim(nameMap.at(StringName::X_DIM), nx);
    netCDF::NcDim yDim = dataGroup.addDim(nameMap.at(StringName::Y_DIM), nx);

    std::vector<netCDF::NcDim> dims = { xDim, yDim };
    for (auto fnNamePair : variableFunctions) {
        const std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims));
        std::vector<double> gathered = gather(data, variableFunctions.at(name));
        var.putVar(gathered.data());
    }
}

void dumpGroup(const std::vector<ElementData>& data, netCDF::NcGroup& headGroup, const NameMap& nameMap)
{
    netCDF::NcGroup metaGroup = headGroup.addGroup(nameMap.at(StringName::METADATA_NODE));
    netCDF::NcGroup dataGroup = headGroup.addGroup(nameMap.at(StringName::DATA_NODE));
    dumpMeta(data, metaGroup, nameMap);
    dumpData(data, dataGroup, nameMap);

}

} /* namespace Nextsim */
