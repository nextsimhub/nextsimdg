/*!
 * @file NetCDFForcings.cpp
 *
 * @date Oct 21, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/NetCDFForcings.hpp"

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

namespace Nextsim {

void NetCDFForcings::setCoordinateIndexFn(const CoordinateIndexFn& fn)
{
    indicesFromLonLat = fn;
}
void NetCDFForcings::setFieldNameLookupFn(const NameLookupFn& fn)
{
    forcingNameFromNSname = fn;
}
void NetCDFForcings::setFileNameLookupFn(const FileLookupFn& fn)
{
    fileNameFn = fn;
}
void NetCDFForcings::setTargetLongitude(const ModelArray& longitudeArray)
{
    targetLon = longitudeArray;
}
void NetCDFForcings::setTargetLatitude(const ModelArray& latitudeArray)
{
    targetLat = latitudeArray;
}

ModelArray NetCDFForcings::getData(const std::string& nsName, const TimePoint& time)
{
    return ModelArray();
}

NetCDFForcings::Buffer NetCDFForcings::getFileIndexData(const std::string& filename, const std::string& fieldName, size_t tIndex)
{
    Buffer data;
    netCDF::NcFile ncFile(filename, netCDF::NcFile::read, netCDF::NcFile::nc4);
    netCDF::NcVar dataVar;
    size_t nx = 0;
    size_t ny = 0;
    // Names of coordinates that should be ignored
    static const std::set<std::string> ignoredFields = {
            "longitude",
            "latitude",
            "time",
            "valid_time",
    };
    size_t nDim;
    for (auto entry : ncFile.getVars()) {
        if (ignoredFields.count(entry.first) == 0 && (fieldName.size() == 0 || entry.first == fieldName)) {
            dataVar = entry.second;
            nDim = dataVar.getDimCount();
            // x and y are assumed to be the last two dimensions
            nx = dataVar.getDim(nDim-1).getSize();
            ny = dataVar.getDim(nDim-2).getSize();
        }
    }
    std::vector<size_t> start(nDim, 0);
    start[0] = tIndex;
    std::vector<size_t> count(nDim, 1);
    count[nDim - 1] = nx;
    count[nDim - 2] = ny;

    data.resize(nx, ny);

    dataVar.getVar(start, count, data.data());
    static const std::string offset_name = "add_offset";
    static const std::string scale_name = "scale_factor";
    auto dataAtts = dataVar.getAtts();
    double a = 1.;
    if (dataAtts.count(scale_name)){
        dataAtts.at(scale_name).getValues(&a);
    }
    double b = 0.;
    if (dataAtts.count(offset_name)) {
        dataAtts.at(offset_name).getValues(&b);
    }
    ncFile.close();

    data *= a;
    data += b;
    return data;
}

NetCDFForcings::Buffer NetCDFForcings::getFileIndexData(const std::string& filename, size_t tIndex)
{
    return getFileIndexData(filename, "", tIndex);
}
ModelArray NetCDFForcings::maFromBuffer(const Buffer& buffer, const ModelArray& fracI, const ModelArray& fracJ)
{
    int nxma = fracI.dimensions()[0];
    int nyma = fracI.dimensions()[1];

    int nxe5 = buffer.rows();
    int nye5 = buffer.cols();
    double ptsPerDegree = nxe5 / 360;
    int nxsrc = nxe5 + 2;
    int nysrc = nye5 + 2;

    Buffer srcBuffer(nxsrc, nysrc);

    srcBuffer(Eigen::seq(1, Eigen::last-1), Eigen::seq(1, Eigen::last-1)) = buffer;
    // Wrap-around columns at the x edges
    srcBuffer(0, Eigen::seq(1, Eigen::last-1)) = buffer(Eigen::last, Eigen::all);
    srcBuffer(Eigen::last, Eigen::seq(1, Eigen::last-1)) = buffer(0, Eigen::all);
    // Duplicate rows at the y edges
    srcBuffer(Eigen::all, 0) = srcBuffer(Eigen::all, 1);
    srcBuffer(Eigen::all, Eigen::last) = srcBuffer(Eigen::all, Eigen::last - 1);

    ModelArray maData(ModelArray::Type::H);
    maData.resize();
    for (size_t j = 0; j < nyma; ++j) {
        for (size_t i = 0; i < nxma; ++i) {
            // Add 1 to account for the halo added above
            double iFloat = fracI(i, j) + 1;
            double jFloat = fracJ(i, j) + 1;
            int ilo = iFloat;
            int ihi = ilo + 1;
            int jlo = jFloat;
            int jhi = jlo + 1;
            double fx = 1 - (iFloat - ilo);
            double fy = 1 - (jFloat - jlo);
            maData(i, j) = fx * fy * srcBuffer(ilo, jlo) +
                    (1 - fx) * fy * srcBuffer(ihi, jlo) +
                    fx * (1 - fy) * srcBuffer(ilo, jhi) +
                    (1 - fx) * (1 - fy) * srcBuffer(ihi, jhi);
        }
    }
    return maData;

}
const NetCDFForcings::Buffer NetCDFForcings::getBufferData(const std::string nsName, const TimePoint& time)
{
    return Buffer();
}

} /* namespace Nextsim */
