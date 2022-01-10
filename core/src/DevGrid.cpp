/*!
 * @file DevGrid.cpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DevGrid.hpp"

#include "/opt/home/include/ncDim.h" // FIXME Remove me
#include "/opt/home/include/ncDouble.h" // FIXME Remove me
#include <cstddef>
#include <ncDim.h>
#include <ncDouble.h>
#include <vector>

// See https://isocpp.org/wiki/faq/pointers-to-members#macro-for-ptr-to-memfn
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

namespace Nextsim {

const std::string DevGrid::ourStructureName = "devgrid";
const std::string xDimName = "x";
const std::string yDimName = "y";
const int DevGrid::nx = 10;
// See https://isocpp.org/wiki/faq/pointers-to-members#array-memfnptrs
// clang-format off
const std::map<std::string, DevGrid::ProgDoubleFn> DevGrid::variableFunctions
    = {    { "hice", &PrognosticData::iceThickness },
           { "cice", &PrognosticData::iceConcentration },
           { "hsnow", &PrognosticData::snowThickness },
           { "sst", &PrognosticData::seaSurfaceTemperature },
           { "sss", &PrognosticData::seaSurfaceSalinity } };
// clang-format on

DevGrid::DevGrid()
    : processedStructureName(ourStructureName)
{
}

DevGrid::~DevGrid() { }

void DevGrid::initMeta(const netCDF::NcGroup& metaGroup) { data.resize(nx * nx); }

void DevGrid::initData(const netCDF::NcGroup& dataGroup)
{
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < nx; ++j) {
            int linearIndex = i * nx + j;
            double hice;
            double cice;
            double hsnow;
            double sst;
            double sss;
            dataGroup.getVar("hice").getVar(std::vector<std::size_t>({ i, j }), &hice);
            dataGroup.getVar("cice").getVar(std::vector<std::size_t>({ i, j }), &cice);
            dataGroup.getVar("hsnow").getVar(std::vector<std::size_t>({ i, j }), &hsnow);
            dataGroup.getVar("sst").getVar(std::vector<std::size_t>({ i, j }), &sst);
            dataGroup.getVar("sss").getVar(std::vector<std::size_t>({ i, j }), &sss);
            // TODO How to store ice temperature data?
            std::array<double, N_ICE_TEMPERATURES> tice = { 0., 0., 0. };
            data = PrognosticData::generate(hice, cice, sst, sss, hsnow, tice);
        }
    }

    for (auto fnNamePair : variableFunctions) {
        std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.getVar(name));
    }
}
void DevGrid::dumpMeta(netCDF::NcGroup& metaGroup) const
{
    // Run the base class output
    IStructure::dumpMeta(metaGroup);
}

void DevGrid::dumpData(netCDF::NcGroup& dataGroup) const
{

    // Create the dimension data, since it has to be in the same group as the
    // data or the parent group
    netCDF::NcDim xDim = dataGroup.addDim(xDimName, nx);
    netCDF::NcDim yDim = dataGroup.addDim(yDimName, nx);

    std::vector<netCDF::NcDim> dims = { xDim, yDim };
    for (auto fnNamePair : variableFunctions) {
        std::string& name = fnNamePair.first;
        netCDF::NcVar var(dataGroup.addVar(name, netCDF::ncDouble, dims));
        std::vector<double> gathered = gather(variableFunctions.at(name));
        var.putVar(gathered.data());
    }
}

std::vector<double> DevGrid::gather(ProgDoubleFn pFunc) const
{
    std::vector<double> gathered(data.size());
    for (int i = 0; i < data.size(); ++i) {
        gathered[i] = CALL_MEMBER_FN(data[i], pFunc);
    }
    return gathered;
}

// Cursor manipulation override functions
int DevGrid::resetCursor()
{
    cursor = data.begin();
    return IStructure::resetCursor();
}
bool DevGrid::validCursor() const { return cursor != data.end(); }
ElementData& DevGrid::cursorData() { return *cursor; }
const ElementData& DevGrid::cursorData() const { return *cursor; }
void DevGrid::incrCursor() { ++cursor; }

} /* namespace Nextsim */
