/*!
 * @author  Einar Olason on 21/07/2026.
 */

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

#include "include/NetCDFUtils.hpp"
#include "include/ParaGridInputs.hpp"

#include "NextsimDynamics.hpp"
#include "include/constants.hpp"
#include "include/indexer.hpp"

namespace Nextsim {

void ParaGridInputs::setData(const TimePoint& time, const std::string& pathSpecIn,
    const std::string& ncLonNameIn, const std::string& ncLatNameIn, const std::string& ncTimeNameIn,
    const std::set<std::string>& forcingsIn,
    const std::set<std::pair<std::string, std::string>>& vectorsIn, const ModelArray& modelLonsIn,
    const ModelArray& modelLatsIn)
{
    currentTime = time;

    pathSpec = pathSpecIn;
    forcings = forcingsIn;
    vectors = vectorsIn;

    ncTimeName = ncTimeNameIn;
    ncLatName = ncLatNameIn;
    ncLonName = ncLonNameIn;

    modelLons = modelLonsIn;
    modelLats = modelLatsIn;

    forcingLonLats = readRawData(currentTime, { ncLatName, ncLonName });

    setWeights();
}

void ParaGridInputs::update(const TimePoint& time)
{
    currentTime = time;

    /* Only actually do something if we find ourselves outside the previous time range. We only
     * check timeRange.after, because time moves forward.
     */
    if (time > timeRange.after) {
        RawDataMap rawDataBefore, rawDataAfter;
        readRawForcing(rawDataBefore, rawDataAfter);

        rotateInputVectors(rawDataBefore);
        rotateInputVectors(rawDataAfter);

        forcingStateBefore = interpolateSpatially(rawDataBefore);
        forcingStateAfter = interpolateSpatially(rawDataAfter);
    }
}

ModelArray ParaGridInputs::getField(const std::string& fieldName)
{
    ModelArray ma;
    ma.reinitialize();

    // Just a linear interpolation between forcing data time steps.
    const FloatType frac = (currentTime - timeRange.before).seconds()
        / (timeRange.after - timeRange.before).seconds();

#pragma omp parallel for
    for (size_t i = 0; i < ma.size(); ++i) {
        ma[i] = frac * forcingStateAfter.data[fieldName][i]
            + (1. - frac) * forcingStateBefore.data[fieldName][i];
    }

    return ma;
}

void ParaGridInputs::setWeights()
{
    // Initialise the size of the weights and indexes
    xi.reinitialize();
    eta.reinitialize();

    ij00.resize(xi.size());
    ij01.resize(xi.size());
    ij10.resize(xi.size());
    ij11.resize(xi.size());

    // Useful aliases
    const auto& lonDimSize = forcingLonLats.dims.at(ncLonName).size();
    const auto& latDimSize = forcingLonLats.dims.at(ncLatName).size();

    // Different methods for Mercator maps and curvilinear grids
    if (lonDimSize == 1 && latDimSize == 1) {
        // Lat and lon are 1D
        setWeights1D();
    } else if (lonDimSize == 2 && latDimSize == 2) {
        // Lat and lon are 2D
        setWeights2D();
    } else {
        throw std::out_of_range(
            "ParaGridInputs::setWeights: Unsupported dimensionality of forcing grid: "
            + std::to_string(lonDimSize) + " & " + std::to_string(latDimSize) + ".\n");
    }
}

void ParaGridInputs::setWeights1D()
{
    // We need to read one of the variables to get the grid ordering
    const auto& choiceVar = readRawData(currentTime, { *forcings.begin() });
    const auto& gridDims = choiceVar.dims.at(*forcings.begin());

    // The axis may be flipped. Probably only the latitude one ... but you never know
    auto forcingLons = forcingLonLats.data[ncLonName];
    auto forcingLats = forcingLonLats.data[ncLatName];

    bool flippedLons = false;
    bool flippedLats = false;
    if (*forcingLons.begin() > *forcingLons.end()) {
        flippedLons = true;
        std::reverse(forcingLons.begin(), forcingLons.end());
    }
    if (*forcingLats.begin() > *forcingLats.end()) {
        flippedLats = true;
        std::reverse(forcingLats.begin(), forcingLats.end());
    }

    // Careful to wrap the longitudes of the model the same as that of the data
    FloatType lon0;
    if (*std::min_element(forcingLons.begin(), forcingLons.end()) < 0.)
        lon0 = 180.;
    else
        lon0 = 360.;

#pragma omp parallel for
    for (size_t i = 0; i < modelLons.size(); ++i) {
        /* Calculate the weights on a local tangent plane, with
         * x = R \cos(\phi) \lambda
         * y = R \phi
         * but R and \cos(\phi) cancel out */

        // This is the target
        const FloatType x = wrapLon(modelLons[i], lon0);
        const FloatType y = modelLats[i];

        // Find the bounding box
        size_t bLon
            = std::upper_bound(forcingLons.begin(), forcingLons.end(), x) - forcingLons.begin();
        size_t bLat
            = std::upper_bound(forcingLats.begin(), forcingLats.end(), y) - forcingLats.begin();

        size_t aLon = bLon - 1;
        size_t aLat = bLat - 1;

        /* Use modulo for periodic boundaries in longitude. Different formulations because aLon
         * maybe negative, but bLon is always positive.
         */
        aLon = (aLon % forcingLons.size() + forcingLons.size()) % forcingLons.size();
        bLon %= forcingLons.size();

        // Now for the weights
        const FloatType x1 = forcingLons[aLon];
        const FloatType x2 = forcingLons[bLon];
        const FloatType y1 = forcingLats[aLat];
        const FloatType y2 = forcingLats[bLat];

        xi[i] = wrapLon(x - x1, lon0) / wrapLon(x2 - x1, lon0);
        eta[i] = (y - y1) / (y2 - y1);

        // If we flipped the axis, we need to point to the right elements in the original vector
        if (flippedLons) {
            aLon = (forcingLons.size() - 1) - aLon;
            bLon = (forcingLons.size() - 1) - bLon;
        }
        if (flippedLats) {
            aLat = (forcingLats.size() - 1) - aLat;
            bLat = (forcingLats.size() - 1) - bLat;
        }

        // Record the bounds
        ij00[i] = indexer(gridDims, { aLon, aLat });
        ij10[i] = indexer(gridDims, { bLon, aLat });
        ij01[i] = indexer(gridDims, { aLon, bLat });
        ij11[i] = indexer(gridDims, { bLon, bLat });
    }
}

void ParaGridInputs::setWeights2D()
{
    /* Just call recursiveBisectSearch for every model grid point. Start off with the full grid. The
     * recursive routine returns false if it didn't find anything at the given level. So a false
     * here means the model point is not in the forcing data set grid.
     */
#pragma omp parallel for
    for (size_t k = 0; k < modelLons.size(); ++k) {
        // Careful with the C-style indexing!
        if (const auto& dims = forcingLonLats.dims.at(ncLonName);
            !recursiveBisectSearch(k, modelLons[k], modelLats[k], 0, dims[0] - 1, 0, dims[1] - 1))
            throw std::out_of_range("ParaGridInputs::setWeights2D: Couldn't find "
                + std::to_string(modelLons[k]) + ", " + std::to_string(modelLats[k])
                + " in the forcing grid.\n");
    }
}

bool ParaGridInputs::recursiveBisectSearch(const size_t k, const FloatType targetLon,
    const FloatType targetLat, const size_t i, const size_t ii, const size_t j, const size_t jj)
{
    // If one of the dimensions has collapsed, then the point is not here(!)
    if (i == ii || j == jj)
        return false;

    // Useful aliases
    const auto& forcingLons = forcingLonLats.data.at(ncLonName);
    const auto& forcingLats = forcingLonLats.data.at(ncLatName);
    const auto& dims = forcingLonLats.dims.at(ncLonName);

    /* Project the corner points onto an orthographic projection, centred on the target. We do the
     * rest of the work in {x,y} coordinates. The target {x,y} is now always at the origin.
     */
    FloatType x00, y00, x10, y10, x01, y01, x11, y11;
    orthographicProjection(forcingLons[indexer(dims, { i, j })],
        forcingLats[indexer(dims, { i, j })], targetLon, targetLat, x00, y00);
    orthographicProjection(forcingLons[indexer(dims, { ii, j })],
        forcingLats[indexer(dims, { ii, j })], targetLon, targetLat, x10, y10);
    orthographicProjection(forcingLons[indexer(dims, { i, jj })],
        forcingLats[indexer(dims, { i, jj })], targetLon, targetLat, x01, y01);
    orthographicProjection(forcingLons[indexer(dims, { ii, jj })],
        forcingLats[indexer(dims, { ii, jj })], targetLon, targetLat, x11, y11);

    // If we're not inside the bounding box, then there's no point in going further
    if (!pointInBoundingBox({ x00, x10, x01, x11 }, { y00, y10, y01, y11 }))
        return false;

    // If the size of the box is one, then we can try to find the local coordinates
    if (ii == i + 1 && jj == j + 1) {
        // Save the index
        ij00[k] = indexer(dims, { i, j });
        ij10[k] = indexer(dims, { ii, j });
        ij01[k] = indexer(dims, { i, jj });
        ij11[k] = indexer(dims, { ii, jj });

        // Try to find local coordinates
        return findLocalCoordinates(k, x00, y00, x10, y10, x01, y01, x11, y11);
    }

    // Bisect and call self to continue searching
    const size_t iHalf = (i + ii) / 2;
    const size_t jHalf = (j + jj) / 2;

    // Search quadrant 00
    if (recursiveBisectSearch(k, targetLon, targetLat, i, iHalf, j, jHalf))
        return true;

    // Search quadrant 10
    if (recursiveBisectSearch(k, targetLon, targetLat, iHalf, ii, j, jHalf))
        return true;

    // Search quadrant 01
    if (recursiveBisectSearch(k, targetLon, targetLat, i, iHalf, jHalf, jj))
        return true;

    // Search quadrant 11
    if (recursiveBisectSearch(k, targetLon, targetLat, iHalf, ii, jHalf, jj))
        return true;

    // Nothing found
    return false;
}

void ParaGridInputs::orthographicProjection(const FloatType lon, const FloatType lat,
    const FloatType lon0, const FloatType lat0, FloatType& x, FloatType& y) const
{
    /* Most of these are used twice, but not all. But anyway, it's easier to read like this, and the
     * compiler should optimise the excessive assignments out, right?
     */
    const FloatType cosPhi = std::cos(lat * PhysicalConstants::deg2rad);
    const FloatType cosPhi0 = std::cos(lat0 * PhysicalConstants::deg2rad);
    const FloatType sinPhi = std::sin(lat * PhysicalConstants::deg2rad);
    const FloatType sinPhi0 = std::sin(lat0 * PhysicalConstants::deg2rad);
    const FloatType cosDeltaLambda = std::cos((lon - lon0) * PhysicalConstants::deg2rad);
    const FloatType sinDeltaLambda = std::sin((lon - lon0) * PhysicalConstants::deg2rad);

    // Projected coordinates
    x = cosPhi * sinDeltaLambda;
    y = cosPhi0 * sinPhi - sinPhi0 * cosPhi * cosDeltaLambda;

    // If the point is outside the projected area, then we move it to the map edge
    if (const FloatType c = sinPhi0 * sinPhi - cosPhi0 * cosPhi * cosDeltaLambda;
        std::cos(c) < 0.) {
        x = std::copysign(1., x);
        y = std::copysign(1., y);
    }
}

bool ParaGridInputs::pointInBoundingBox(
    const std::vector<FloatType>& xCorners, const std::vector<FloatType>& yCorners)
{
    const FloatType xMin = *std::min_element(xCorners.begin(), xCorners.end());
    const FloatType xMax = *std::max_element(xCorners.begin(), xCorners.end());

    const FloatType yMin = *std::min_element(yCorners.begin(), yCorners.end());
    const FloatType yMax = *std::max_element(yCorners.begin(), yCorners.end());

    // The target point is at the origin
    return 0. >= xMin && 0. <= xMax && 0. >= yMin && 0. <= yMax;
}

bool ParaGridInputs::findLocalCoordinates(const size_t k, const FloatType x00, const FloatType y00,
    const FloatType x10, const FloatType y10, const FloatType x01, const FloatType y01,
    const FloatType x11, const FloatType y11)
{
    /*
     * Cell corners:
     *
     *       p01 -------- p11
     *        |            |
     *        |            |
     *       p00 -------- p10
     *
     * Local coordinates:
     *
     *       eta = 1
     *          ^
     *          |
     *          |
     *       eta = 0
     *
     *       xi = 0 ---> xi = 1
     */

    // Initial guess.
    // The center of the cell is usually a reasonable starting point.
    xi[k] = 0.5;
    eta[k] = 0.5;

    /* Newton iteration:
     *     F(xi, eta) = mapping(xi, eta) - query_point
     */
    for (int iteration = 0; iteration < 20; ++iteration) {
        // What is the right tolerance for both single and double? I guess this is fine.
        constexpr FloatType tolerance = 10 * std::numeric_limits<FloatType>::epsilon();

        // Bilinear mapping from local coordinates to physical coordinates.
        const FloatType N00 = (1.0 - xi[k]) * (1.0 - eta[k]);
        const FloatType N10 = xi[k] * (1.0 - eta[k]);
        const FloatType N01 = (1.0 - xi[k]) * eta[k];
        const FloatType N11 = xi[k] * eta[k];

        const FloatType xp = N00 * x00 + N10 * x10 + N01 * x01 + N11 * x11;
        const FloatType yp = N00 * y00 + N10 * y10 + N01 * y01 + N11 * y11;

        // Convergence test - target is at origin
        if (std::sqrt(xp * xp + yp * yp) < tolerance)
            break;

        /* Jacobian:
         *     [ dx/dxi   dx/deta ]
         * J = [                  ]
         *     [ dy/dxi   dy/deta ]
         */
        const FloatType dx_dxi = (1.0 - eta[k]) * (x10 - x00) + eta[k] * (x11 - x01);
        const FloatType dx_deta = (1.0 - xi[k]) * (x01 - x00) + xi[k] * (x11 - x10);
        const FloatType dy_dxi = (1.0 - eta[k]) * (y10 - y00) + eta[k] * (y11 - y01);
        const FloatType dy_deta = (1.0 - xi[k]) * (y01 - y00) + xi[k] * (y11 - y10);

        const FloatType determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;

        // Check if the cell is degenerate
        if (std::abs(determinant) < tolerance)
            return false;

        /* Solve:
         *     J [d_xi ] = -F
         *       [d_eta]
         */
        xi[k] -= (xp * dy_deta - dx_deta * yp) / determinant;
        eta[k] += (dy_dxi * xp - dx_dxi * yp) / determinant;

        // If the solution is far outside the cell, then this is probably not the correct cell.
        if (xi[k] < -0.1 || xi[k] > 1.1 || eta[k] < -0.1 || eta[k] > 1.1)
            return false;
    }

    // Check whether the point is actually inside the cell.
    if (!(xi[k] >= 0. && xi[k] <= 1. && eta[k] >= 0. && eta[k] <= 1.))
        return false;

    return true;
}

ModelState ParaGridInputs::interpolateSpatially(const RawDataMap& rawData)
{
    ModelState state;
    for (const auto& dataPair : rawData.data) {
        // Structured bindings and omp don't mesh
        const std::string& name = dataPair.first;
        const std::vector<FloatType>& data = dataPair.second;

        state.data[name].reinitialize();
#pragma omp parallel for
        for (size_t i = 0; i < state.data.at(name).size(); ++i) {
            const FloatType f00 = data[ij00[i]];
            const FloatType f10 = data[ij10[i]];
            const FloatType f01 = data[ij01[i]];
            const FloatType f11 = data[ij11[i]];

            const FloatType N00 = (1.0 - xi[i]) * (1.0 - eta[i]);
            const FloatType N10 = xi[i] * (1.0 - eta[i]);
            const FloatType N01 = (1.0 - xi[i]) * eta[i];
            const FloatType N11 = xi[i] * eta[i];

            state.data.at(name)[i] = N00 * f00 + N10 * f10 + N01 * f01 + N11 * f11;
        }
    }

    return state;
}

void ParaGridInputs::rotateInputVectors(RawDataMap& rawData)
{
    const auto originalRotation = rawData.data;
    for (const auto& [first, second] : vectors) {
        vectorRotationLogic(originalRotation.at(first), originalRotation.at(second),
            rawData.data.at(first), rawData.data.at(second));
    }
}

void ParaGridInputs::vectorRotationLogic(const std::vector<FloatType>& vectorIn1st,
    const std::vector<FloatType>& vectorIn2nd, std::vector<FloatType>& vectorOut1st,
    std::vector<FloatType>& vectorOut2nd)
{
    // Nothing done yet
    vectorOut1st = vectorIn1st;
    vectorOut2nd = vectorIn2nd;
}

void ParaGridInputs::readRawForcing(RawDataMap& rawDataBefore, RawDataMap& rawDataAfter)
{
    /* First we find the correct time, time slice, and files. For multi-file datasets we just pick
     * any (first) file/variable.
     */
    const std::string fileName = formatFileName(currentTime, *forcings.begin());
    size_t targetTIndexAfter, targetTIndexBefore;
    try {
        netCDF::NcFile ncFile(fileName, netCDF::NcFile::read);

        // Read the time axis
        netCDF::NcDim timeDim = ncFile.getDim(ncTimeName);
        // Read the time variable
        netCDF::NcVar timeVar = ncFile.getVar(ncTimeName);
        // Get the time axis as a vector. We use double here, because Duration expects double.
        std::vector<double> timeVec(timeDim.getSize());
        timeVar.getVar(timeVec.data());

        // Time units nonsense
        std::string unitStr, timeUnit, sinceKeyWord, timeOrigin;
        timeVar.getAtt("units").getValues(unitStr);
        std::stringstream ss(unitStr);
        ss >> timeUnit >> sinceKeyWord >> timeOrigin;

        // Get a TimePoint for the origin
        const auto timePointOrigin = TimePoint(timeOrigin);

        // Multiply the time axis to get seconds
        double multiplier;
        if (timeUnit == "seconds")
            multiplier = 1.;
        else if (timeUnit == "minutes")
            multiplier = 60.;
        else if (timeUnit == "hour" || timeUnit == "hours")
            multiplier = 3600.;
        else if (timeUnit == "days")
            multiplier = 24. * 3600.;
        else
            throw std::runtime_error("ParaGridInputs::readRawForcing(): unsupported time unit "
                + timeUnit + " in '" + unitStr + "'.\n");

        std::transform(timeVec.begin(), timeVec.end(), timeVec.begin(),
            [multiplier](const double t) { return t * multiplier; });

        // Because currentTime can't be captured in the lambda function
        const TimePoint& time = currentTime;

        // Get the index of the first TimePoint greater than the target.
        targetTIndexAfter = std::find_if(timeVec.begin(), timeVec.end(),
                                [time, timePointOrigin](const double t) {
                                    return TimePoint(timePointOrigin, Duration(t)) > time;
                                })
            - timeVec.begin();

        targetTIndexBefore = targetTIndexAfter - 1;
        timeRange.before = TimePoint(timePointOrigin, Duration(timeVec[targetTIndexBefore]));
        timeRange.after = TimePoint(timePointOrigin, Duration(timeVec[targetTIndexAfter]));

        /* We need to check if targetTIndexAfter is actually pointing at the right time, or if we
         * need to go on to the next file. Lucky for us, std::find_if returns timeVec.end() if it
         * finds nothing, so targetTIndexBefore is always right.
         */
        if (targetTIndexAfter == timeVec.size()) {
            // Assume time is increasing at a constant rate(!)
            timeRange.after = timeRange.before + Duration(timeVec[1] - timeVec[0]);
            targetTIndexAfter = 0;
        }

        // Sanity check. Not really needed.
        if (targetTIndexAfter < 0 || targetTIndexBefore < 0 || targetTIndexAfter >= timeVec.size()
            || targetTIndexBefore >= timeVec.size())
            throw std::out_of_range(
                "ParaGridInputs::readRawForcing::Target time index is out of range "
                "- how could this happen?\n");

        ncFile.close();
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + fileName;
        throw std::runtime_error(ncWhat);
    }

    // Read the data
    rawDataBefore = readRawData(timeRange.before, forcings, targetTIndexBefore);
    rawDataAfter = readRawData(timeRange.after, forcings, targetTIndexAfter);
}

ParaGridInputs::RawDataMap ParaGridInputs::readRawData(
    const TimePoint& time, const std::set<std::string>& fields, const size_t timeIndex) const
{
    RawDataMap data;
    std::string fileName;

    // Loop over the variable names at the top, because we may have one variable per file.
    for (const std::string& varName : fields) {
        try {

            // We may need to be more careful selecting the right file if we're reading lat/lon info
            if (varName == ncLonName || varName == ncLatName)
                fileName = formatFileName(time, *forcings.begin());
            else
                fileName = formatFileName(time, varName);

            netCDF::NcFile ncFile(fileName, netCDF::NcFile::read);

            // Don't try to read non-existent data
            if (ncFile.getVars().count(varName) == 0)
                continue;

            netCDF::NcVar var = ncFile.getVar(varName);

            std::vector<netCDF::NcDim> dims = var.getDims();
            std::vector<size_t> start(dims.size(), 0);
            std::vector<size_t> count(dims.size());
            for (int i = 0; i < count.size(); ++i)
                count[i] = dims[i].getSize();

            // Pick the time slice if we have a time axis
            for (int i = 0; i < dims.size(); ++i) {
                if (dims[i].getName() == ncFile.getVar(ncTimeName).getDims()[0].getName()) {
                    start[i] = timeIndex;
                    count[i] = 1;
                    continue;
                }
                // Push all non-time dimension counts to the output data map
                data.dims[varName].push_back(count[i]);
            }
            // Needs to be reversed because of netCDF shenanigans
            std::reverse(data.dims[varName].begin(), data.dims[varName].end());

            // Resize and read!
            data.data[varName] = std::vector<FloatType>(std::accumulate(
                count.begin(), count.end(), static_cast<size_t>(1), std::multiplies<>()));
            var.getVar(start, count, data.data[varName].data());

        } catch (const netCDF::exceptions::NcException& nce) {
            std::string ncWhat(nce.what());
            ncWhat += ": " + fileName;
            throw std::runtime_error(ncWhat);
        }
    }
    return data;
}

}