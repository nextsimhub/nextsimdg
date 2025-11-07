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
    double missing = 0.;
    return getFileIndexData(filename, fieldName, tIndex, missing);
}

NetCDFForcings::Buffer NetCDFForcings::getFileIndexData(const std::string& filename, const std::string& fieldName, size_t tIndex, double& missing)
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
    static const std::string offsetName = "add_offset";
    static const std::string scaleName = "scale_factor";
    static const std::string missingName = "missing_value";
    static const std::string fillValueName = "_FillValue";
    auto dataAtts = dataVar.getAtts();
    double a = 1.;
    if (dataAtts.count(scaleName)){
        dataAtts.at(scaleName).getValues(&a);
    }
    double b = 0.;
    if (dataAtts.count(offsetName)) {
        dataAtts.at(offsetName).getValues(&b);
    }
    missing = 0.;
    bool hasMissing = false;
    if (dataAtts.count(missingName)) {
        dataAtts.at(missingName).getValues(&missing);
        hasMissing = true;
    } else if (dataAtts.count(fillValueName)) {
        dataAtts.at(fillValueName).getValues(&missing);
        hasMissing = true;
    }
    ncFile.close();

    data *= a;
    data += b;

    // Don't scale and offset the missing value if one was not set.
    if (hasMissing) {
        missing *= a;
        missing += b;
    }
    return data;
}

NetCDFForcings::Buffer NetCDFForcings::getFileIndexData(const std::string& filename, size_t tIndex)
{
    return getFileIndexData(filename, "", tIndex);
}

NetCDFForcings::Buffer NetCDFForcings::getFileIndexData(const std::string& filename, size_t tIndex, double& missing)
{
    return getFileIndexData(filename, "", tIndex, missing);
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

/*!
 * @brief Interpolates from the native forcing grid to the model grid.
 *
 * @details Interpolates from the grid on which the forcing value is natively
 * stored to the model grid. Missing data is ignored in the interpolation, and
 * a missing data value of 0. is taken to mean that there is no missing data
 * value provided, as there are few circumstances where a value of 0. would not
 * be a possible valid value. If the original data does use 0. as a missing
 * data value, then the data should be re-masked before this function is
 * called.
 *
 * @param buffer The forcing data on its native grid
 * @param iFrac The x index (including fractional part) on the forcing grid for
 *     each point on the target model grid.
 * @param jFrac The y index (including fractional part) on the forcing grid for
 *     each point on the target model grid.
 * @param missing The value indicating missing data. A value of 0 is assumed to mean there is no missing data value
 */
ModelArray NetCDFForcings::maFromBuffer(const Buffer& buffer, const ModelArray& fracI, const ModelArray& fracJ, double missing)
{
    const size_t nSamples = 4;
    const std::array<int, nSamples> xOffsets = { 0, 1, 0, 1 };
    const std::array<int, nSamples> yOffests = { 0, 0, 1, 1 };
    int nxma = fracI.dimensions()[0];
    int nyma = fracI.dimensions()[1];

    int nxbu = buffer.rows();
    int nybu = buffer.cols();
    ModelArray maData(ModelArray::Type::H);
    maData.resize();
    for (size_t j = 0; j < nyma; ++j) {
        for (size_t i = 0; i < nxma; ++i) {
            double iFloat = fracI(i, j);
            double jFloat = fracJ(i, j);
            // Separate out the integral and fractional parts
            int ilo = iFloat;
            int jlo = jFloat;
            double fx = 1 - (iFloat - ilo);
            double fy = 1 - (jFloat - jlo);
            const std::array<double, nSamples> weights = { fx*fy, (1-fx)*fy, fx*(1-fy), (1-fx)*(1-fy) };
            double vAcc = 0.;
            double wAcc = 0.;
            for (size_t s = 0; s < nSamples; ++s) {
                if (ilo >= 0 && ilo+1 < nxbu && jlo >= 0 && jlo+1 < nybu) {
                    double v = buffer(ilo, jlo);
                    if (missing == 0 || v != missing) {
                        vAcc += v * weights[s];
                        wAcc += weights[s];
                    }
                }
                maData(i, j) = (wAcc != 0.) ? vAcc / wAcc : missing;
            }
        }
    }

    return maData;
}
const NetCDFForcings::Buffer NetCDFForcings::getBufferData(const std::string nsName, const TimePoint& time)
{
    return Buffer();
}

} /* namespace Nextsim */
