/*!
 * @author  Einar Olason on 21/07/2026.
 */

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

#include "include/NetCDFUtils.hpp"
#include "include/ParaGridInputs.hpp"
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

    /* Here we:
     *  - Read the dimensions in the netCDF file
     *  - Read the lon/lat coordinates
     *  - Calculate the weights
     *  - Calculate a reduced domain
     *  - Re-read the lon/lat coordinates
     * Now, everything should just work on the reduced coordinates. This reduces memory usage, and
     * is also good for MPI).
     */
    readDims();
    forcingLonLats = readRawData<double>(currentTime, { ncLatName, ncLonName });
    setWeights();
    tightenGrid();
    forcingLonLats = readRawData<double>(currentTime, { ncLatName, ncLonName });

    // Different methods for Mercator maps and curvilinear grids
    VectorRotator::orientation orient;
    if (lonLat1D)
        orient = VectorRotator::orientation::EAST_NORTH;
    else
        orient = VectorRotator::orientation::GRID;

    rotator = std::make_unique<VectorRotator>(
        gridDims, forcingLonLats.at(ncLonName), forcingLonLats.at(ncLatName), orient);
}

void ParaGridInputs::tightenGrid()
{
    gridStart = { std::numeric_limits<size_t>::max(), std::numeric_limits<size_t>::max() };
    std::vector<size_t> gridEnd = { 0, 0 };

    // Loop over all the points in the corner lists to find the grid boundaries
    for (const auto* cornerPtr : { &ij00, &ij01, &ij10, &ij11 }) {
#pragma omp parallel for
        for (const auto& point : *cornerPtr) {
            const auto ij = deIndexer(gridDims, point);
            for (size_t k = 0; k < ij.size(); k++) {
                gridStart[k] = std::min(gridStart[k], ij[k]);
                gridEnd[k] = std::max(gridEnd[k], ij[k]);
            }
        }
    }

    // Update grid dimensions, now that we have start and end values
    // Careful with one-off!
    const std::vector<size_t> oldDims = gridDims;
    gridDims[0] = gridEnd[0] - gridStart[0] + 1;
    gridDims[1] = gridEnd[1] - gridStart[1] + 1;

    // Loop again over the corner lists to shift the coordinates
    for (auto* cornerPtr : { &ij00, &ij01, &ij10, &ij11 }) {
#pragma omp parallel for
        for (auto& point : *cornerPtr) {
            const auto ij = deIndexer(oldDims, point);
            point = indexer(gridDims, { ij[0] - gridStart[0], ij[1] - gridStart[1] });
        }
    }
}

void ParaGridInputs::readDims()
{
    // Get the grid dimensions from the frist variable in forcings
    const std::string& varName = *forcings.begin();
    const std::string fileName = formatFileName(currentTime, varName);
    try {
        const netCDF::NcFile ncFile(fileName, netCDF::NcFile::read);

        const std::vector<netCDF::NcDim> dims = ncFile.getVar(varName).getDims();
        const std::string timeDimName = ncFile.getVar(ncTimeName).getDims()[0].getName();
        for (const auto& dim : dims) {
            if (dim.getName() != timeDimName)
                gridDims.push_back(dim.getSize());
        }

        // Needs to be reversed because of netCDF shenanigans
        std::reverse(gridDims.begin(), gridDims.end());
        gridStart.assign(gridDims.size(), 0);

        // Read the dimensions of lat and long variables
        const std::vector<netCDF::NcDim> lonDims = ncFile.getVar(ncLonName).getDims();
        const std::vector<netCDF::NcDim> latDims = ncFile.getVar(ncLatName).getDims();

        if (latDims.size() == 1 && lonDims.size() == 1) {
            lonLat1D = true;
            if (lonDims[0].getSize() != gridDims[0] || latDims[0].getSize() != gridDims[1])
                throw std::runtime_error(
                    "ParaGridInputs::readDims: Inconsistent dimension sizes for " + varName
                    + " and longitude and latitude variables: [" + std::to_string(gridDims[0]) + ","
                    + std::to_string(gridDims[1]) + "] and [" + std::to_string(lonDims[0].getSize())
                    + "," + std::to_string(latDims[0].getSize()) + "] respectively.\n");
        } else if (latDims.size() == 2 && lonDims.size() == 2) {
            lonLat1D = false;
            if (latDims[1].getSize() != gridDims[0] || latDims[0].getSize() != gridDims[1])
                throw std::runtime_error(
                    "ParaGridInputs::readDims: Inconsistent dimension sizes for " + varName
                    + " and longitude and latitude variables: [" + std::to_string(gridDims[0]) + ","
                    + std::to_string(gridDims[1]) + "] and [" + std::to_string(latDims[0].getSize())
                    + "," + std::to_string(latDims[1].getSize()) + "] respectively.\n");
            if (lonDims[1].getSize() != gridDims[0] || lonDims[0].getSize() != gridDims[1])
                throw std::runtime_error(
                    "ParaGridInputs::readDims: Inconsistent dimension sizes for " + varName
                    + " and longitude and latitude variables: [" + std::to_string(gridDims[0]) + ","
                    + std::to_string(gridDims[1]) + "] and [" + std::to_string(lonDims[0].getSize())
                    + "," + std::to_string(lonDims[1].getSize()) + "] respectively.\n");
        } else {
            throw std::runtime_error("ParaGridInputs::readDims: Inconsistent dimension size for "
                + ncLonName + " and " + ncLatName + " " + std::to_string(lonDims.size()) + " and "
                + std::to_string(latDims.size()) + " respectively.\n");
        }
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + fileName;
        throw std::runtime_error(ncWhat);
    }
}

void ParaGridInputs::update(const TimePoint& time)
{
    currentTime = time;

    /* Only actually do something if we find ourselves outside the previous time range. We only
     * check timeRange.after, because time moves forward.
     */
    if (time > timeRange.after) {
        RawDataMap<FloatType> rawDataBefore, rawDataAfter;
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

    // Different methods for Mercator maps and curvilinear grids
    if (lonLat1D)
        setWeights1D();
    else
        setWeights2D();
}

void ParaGridInputs::setWeights1D()
{
    // Useful alias
    auto& forcingLons = forcingLonLats.at(ncLonName);

    // The latitude axis may be flipped, so it can't be a reference.
    auto forcingLats = forcingLonLats.at(ncLatName);

    bool flippedLats = false;
    if (*forcingLats.begin() > *forcingLats.end()) {
        flippedLats = true;
        std::reverse(forcingLats.begin(), forcingLats.end());
    }

    // Careful to wrap the longitudes of the model the same as that of the data
    const FloatType lon0 = *std::min_element(forcingLons.begin(), forcingLons.end());

#pragma omp parallel for
    for (size_t i = 0; i < modelLons.size(); ++i) {
        /* Calculate the weights on a local tangent plane, with
         * x = R \cos(\phi) \lambda
         * y = R \phi
         * but R and \cos(\phi) cancel out */

        // This is the target. It can be FloatType, because we don't use any trigometric functions
        const FloatType x = wrapLon(modelLons[i], lon0);
        const FloatType y = modelLats[i];

        // Find the bounding box
        size_t bLon
            = std::upper_bound(forcingLons.begin(), forcingLons.end(), x) - forcingLons.begin();
        size_t bLat
            = std::upper_bound(forcingLats.begin(), forcingLats.end(), y) - forcingLats.begin();

        // Just a quick check on latitude.
        // TODO: Check also longitude, taking periodicity into account.
        if (bLat == forcingLats.size())
            throw std::out_of_range("ParaGridInputs::setWeights1D: Couldn't find "
                + std::to_string(x) + ", " + std::to_string(y) + " in the forcing grid.\n");

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
        if (!recursiveBisectSearch(
                k, modelLons[k], modelLats[k], 0, gridDims[0] - 1, 0, gridDims[1] - 1))
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
    const auto& forcingLons = forcingLonLats.at(ncLonName);
    const auto& forcingLats = forcingLonLats.at(ncLatName);

    /* Project the corner points onto an orthographic projection, centred on the target. We do the
     * rest of the work in {x,y} coordinates. The target {x,y} is now always at the origin.
     */
    FloatType x00, y00, x10, y10, x01, y01, x11, y11;
    orthographicProjection(forcingLons[indexer(gridDims, { i, j })],
        forcingLats[indexer(gridDims, { i, j })], targetLon, targetLat, x00, y00);
    orthographicProjection(forcingLons[indexer(gridDims, { ii, j })],
        forcingLats[indexer(gridDims, { ii, j })], targetLon, targetLat, x10, y10);
    orthographicProjection(forcingLons[indexer(gridDims, { i, jj })],
        forcingLats[indexer(gridDims, { i, jj })], targetLon, targetLat, x01, y01);
    orthographicProjection(forcingLons[indexer(gridDims, { ii, jj })],
        forcingLats[indexer(gridDims, { ii, jj })], targetLon, targetLat, x11, y11);

    // If we're not inside the bounding box, then there's no point in going further
    if (!pointInBoundingBox({ x00, x10, x01, x11 }, { y00, y10, y01, y11 }))
        return false;

    // If the size of the box is one, then we can try to find the local coordinates
    if (ii == i + 1 && jj == j + 1) {
        // Save the index
        ij00[k] = indexer(gridDims, { i, j });
        ij10[k] = indexer(gridDims, { ii, j });
        ij01[k] = indexer(gridDims, { i, jj });
        ij11[k] = indexer(gridDims, { ii, jj });

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

void ParaGridInputs::orthographicProjection(const double lon, const double lat,
    const FloatType lon0, const FloatType lat0, FloatType& x, FloatType& y)
{
    /* Most of these are used twice, but not all. But anyway, it's easier to read like this, and the
     * compiler should optimise the excessive assignments out, right?
     */
    const double cosPhi = std::cos(radians(lat));
    const double cosPhi0 = std::cos(radians(lat0));
    const double sinPhi = std::sin(radians(lat));
    const double sinPhi0 = std::sin(radians(lat0));
    const double cosDeltaLambda = std::cos(radians(lon - lon0));
    const double sinDeltaLambda = std::sin(radians(lon - lon0));

    // Projected coordinates
    x = cosPhi * sinDeltaLambda;
    y = cosPhi0 * sinPhi - sinPhi0 * cosPhi * cosDeltaLambda;

    // If the point is outside the projected area, then we move it to the map edge
    if (const double c = sinPhi0 * sinPhi - cosPhi0 * cosPhi * cosDeltaLambda; std::cos(c) < 0.) {
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

ModelState ParaGridInputs::interpolateSpatially(const RawDataMap<FloatType>& rawData)
{
    ModelState state;
    for (const auto& dataPair : rawData) {
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

void ParaGridInputs::rotateInputVectors(RawDataMap<FloatType>& rawData)
{
    // Useful alias
    const auto& forcingLats = forcingLonLats.at(ncLatName);

    for (const auto& [uName, vName] : vectors) {
        auto& uData = rawData.at(uName);
        auto& vData = rawData.at(vName);

        rotator->fromParametricMesh(uData, vData);

        /* We may have a lat/lon dataset with vector data at the pole (ERA5)!
         * Now that we've rotated the vectors, the proper value at the pole can reasonably be
         * interpolated as the mean of all surrounding values.
         * The longitude limit of 89.9 degrees corresponds to sin(89.9) = 0.999998 (five nines), so
         * if we get closer to the pole thant his, then we'll have problems with the vector rotator
         * (assuming double precision).
         */
        if (lonLat1D && *std::max_element(forcingLats.begin(), forcingLats.end()) >= 89.9) {
            rotator->fixLonLatPole(uData, vData, forcingLats);
        }
    }
}

void ParaGridInputs::readRawForcing(
    RawDataMap<FloatType>& rawDataBefore, RawDataMap<FloatType>& rawDataAfter)
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
    rawDataBefore = readRawData<FloatType>(timeRange.before, forcings, targetTIndexBefore);
    rawDataAfter = readRawData<FloatType>(timeRange.after, forcings, targetTIndexAfter);
}

template <typename T>
ParaGridInputs::RawDataMap<T> ParaGridInputs::readRawData(
    const TimePoint& time, const std::set<std::string>& fields, const size_t timeIndex) const
{
    RawDataMap<T> data;
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

            /* Populate start and count, based on the already established gridStart and gridDims,
             * while taking the time dimension into account.
             * NB! j needs to run backwards because of netCDF shenanigans
             * NB! If we're reading lon/lat from a Mercator map we can assume lon is the first and
             * lat the second dimension.
             */
            std::vector<size_t> start, count;
            if (lonLat1D && varName == ncLonName) {
                start.push_back(gridStart[0]);
                count.push_back(gridDims[0]);
            } else if (lonLat1D && varName == ncLatName) {
                start.push_back(gridStart[1]);
                count.push_back(gridDims[1]);
            } else {
                size_t j = gridDims.size() - 1;
                for (const auto& dim : dims) {
                    if (dim.getName() == ncFile.getVar(ncTimeName).getDims()[0].getName()) {
                        start.push_back(timeIndex);
                        count.push_back(1);
                    } else {
                        start.push_back(gridStart[j]);
                        count.push_back(gridDims[j]);
                        --j;
                    }
                }
            }

            // Resize and read!
            data[varName] = std::vector<T>(std::accumulate(
                count.begin(), count.end(), static_cast<size_t>(1), std::multiplies<>()));
            readNetCDFVar(var, start, count, data.at(varName).data());

        } catch (const netCDF::exceptions::NcException& nce) {
            std::string ncWhat(nce.what());
            ncWhat += ": " + fileName;
            throw std::runtime_error(ncWhat);
        }
    }
    return data;
}

}