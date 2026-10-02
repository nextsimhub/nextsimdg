/*!
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#ifndef NEXTSIM_DG_PARAGRIDINPUTS_HPP
#define NEXTSIM_DG_PARAGRIDINPUTS_HPP

#include "ModelArray.hpp"
#include "ModelState.hpp"
#include "Time.hpp"

#include <regex>
#include <set>
#include <string>

namespace Nextsim {

class ParaGridInputs {

public:
    ParaGridInputs() = default;

    /*!
     *@brief Initialise object parameters and calculate interpolation weights. To be called by
     * parent's setData function.
     *
     * @param time Initial time
     * @param pathSpecIn Pathspec for the input file(s). This will be formatted with any time
     * formats (e.g. %Y or %m) for the current time. For datasets comprising of multiple files with
     * one file per variable, the key %v will be replaced by the variable name.
     * @param ncLonNameIn Name of the longitude variable in the netCDF file
     * @param ncLatNameIn Name of the latitude variable in the netCDF file
     * @param ncTimeNameIn Name of the time variable in the netCDF file
     * @param forcingsIn Set of forcing variables to read
     * @param vectorsIn Set of vector variables to read
     * @param modelLonsIn Model longitudes
     * @param modelLatsIn Model latitudes
     */
    void setData(const TimePoint& time, const std::string& pathSpecIn,
        const std::string& ncLonNameIn, const std::string& ncLatNameIn,
        const std::string& ncTimeNameIn, const std::set<std::string>& forcingsIn,
        const std::set<std::pair<std::string, std::string>>& vectorsIn,
        const ModelArray& modelLonsIn, const ModelArray& modelLatsIn);

    /*!
     * @brief Update the forcing data for the current time step.
     *
     * @param time Current model time
     */
    void update(const TimePoint& time);

    /*!
     * @brief Get the forcing data for a specific field interpolated in time to currentTime.
     *
     * @param fieldName Name of the field to retrieve
     * @return Forcing data for the specified field
     */
    [[nodiscard]] ModelArray getField(const std::string& fieldName);

    /*!
     * @brief Check if a field can be read in from the forcing data set.
     *
     * @param fieldName Name of the field to check
     * @return True if the field is valid, false otherwise
     */
    [[nodiscard]] bool isValid(const std::string& fieldName) const
    {
        return forcingStateBefore.data.count(fieldName) && forcingStateAfter.data.count(fieldName);
    }

private:
    // Useful structs
    typedef struct {
        std::map<std::string, std::vector<size_t>> dims;
        std::map<std::string, std::vector<FloatType>> data;
    } RawDataMap;

    struct {
        /* The initialisation is important, because init time should always be larger than
         * timeRange.after
         */
        TimePoint before = std::numeric_limits<TimePoint>::min(),
                  after = std::numeric_limits<TimePoint>::min();
    } timeRange;

    // Weights and coordinate pointers for the bi-linear interpolation
    ModelArray xi, eta;
    std::vector<size_t> ij00, ij01, ij10, ij11;

    // Some (hopefully) self explanatory state variables
    TimePoint currentTime;
    ModelState forcingStateBefore, forcingStateAfter;
    std::string pathSpec, ncTimeName, ncLonName, ncLatName;
    ModelArray modelLons, modelLats;
    RawDataMap forcingLonLats;
    std::set<std::string> forcings;
    std::set<std::pair<std::string, std::string>> vectors;

    // Basic weight-setting functions
    void setWeights();
    void setWeights1D();
    void setWeights2D();

    /* Does a recursive search for the grid cell with corners i, ii, j, and jj, in which the point
     * {targetLon, targetLat} is found.
     */
    bool recursiveBisectSearch(size_t k, FloatType targetLon, FloatType targetLat, size_t i,
        size_t ii, size_t j, size_t jj);

    // Project {lon, lat} onto the orthographic coordinates {x,y} with {lon0, lat0} at its centre.
    void orthographicProjection(FloatType lon, FloatType lat, FloatType lon0, FloatType lat0,
        FloatType& x, FloatType& y) const;

    /* Do a quick axis-aligned bounding box check to see if a point is (probably) in the grid cell.
     * This is a necessary condition, not a sufficient one.
     */
    [[nodiscard]] bool pointInBoundingBox(
        const std::vector<FloatType>& xCorners, const std::vector<FloatType>& yCorners);

    /* Find the local coordinates and weights {xi, eta} for a bi-linear interpolation on a
     * curvilinear grid
     */
    [[nodiscard]] bool findLocalCoordinates(size_t k, FloatType x00, FloatType y00, FloatType x10,
        FloatType y10, FloatType x01, FloatType y01, FloatType x11, FloatType y11);

    // Apply the weights to do a bi-linear interpolation
    [[nodiscard]] ModelState interpolateSpatially(const RawDataMap& rawData);

    // Rotate the vectors from the input to model grid
    void rotateInputVectors(RawDataMap& rawData);

    // A placeholder for the actual vector rotation logic
    void vectorRotationLogic(const std::vector<FloatType>& vectorIn1st,
        const std::vector<FloatType>& vectorIn2nd, std::vector<FloatType>& vectorOut1st,
        std::vector<FloatType>& vectorOut2nd);

    // Read the forcing listed in ``forcings`` at times bracketing ``currentTime``.
    void readRawForcing(RawDataMap& rawDataBefore, RawDataMap& rawDataAfter);

    // The function that actually reads data from the netCDF file
    [[nodiscard]] RawDataMap readRawData(
        const TimePoint& time, const std::set<std::string>& fields, size_t timeIndex = 0) const;

    // Wrap longitudes, so that lon0 is the largest value (usually either [-180 180] or [0 360]
    [[nodiscard]] FloatType wrapLon(const FloatType lon, const FloatType lon0 = 0.) const
    {
        // Shift the value relative to the lower bound so the range starts at 0
        const FloatType shifted = lon - lon0;

        // Wrap safely to [0, 360) using the branchless double fmod trick
        const FloatType wrapped = std::fmod(std::fmod(shifted, 360._ft) + 360._ft, 360._ft);

        // Shift back to the original custom baseline
        return wrapped + lon0;
    }

    // Format the pathspec, replacing %v with a given variable name and using TimePoint.format()
    [[nodiscard]] std::string formatFileName(
        const TimePoint& time, const std::string& variable = "NONE") const
    {
        const std::regex pattern("%v");
        return time.format(std::regex_replace(pathSpec, pattern, variable));
    }
};

}

#endif // NEXTSIM_DG_PARAGRIDINPUTS_HPP
