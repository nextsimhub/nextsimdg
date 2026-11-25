/*!
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#ifndef GEOREFERENCEDINPUT_HPP
#define GEOREFERENCEDINPUT_HPP

#include "include/IForcingInput.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"
#include "include/VectorRotator.hpp"

#include <ncFile.h>
#include <regex>
#include <set>
#include <string>

namespace Nextsim {

class GeoreferencedInput : public IForcingInput {

public:
    GeoreferencedInput() = default;

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
    void initialise(const TimePoint& time, const std::string& pathSpecIn,
        const std::string& ncLonNameIn, const std::string& ncLatNameIn,
        const std::string& ncTimeNameIn, const std::set<std::string>& forcingsIn,
        const std::set<std::pair<std::string, std::string>>& vectorsIn,
        const ModelArray& modelLonsIn, const ModelArray& modelLatsIn) override;

    /*!
     * @brief Update the forcing data for the current time step.
     *
     * @param time Current model time
     */
    void update(const TimePoint& time) override;

    /*!
     * @brief Get the forcing data for a specific field interpolated in time to currentTime.
     *
     * @param fieldName Name of the field to retrieve
     * @param destination The array to populate with the forcing data.
     */
    void getField(const std::string& fieldName, ModelArrayAuto& destination) const override;

    /*!
     * @brief Check if a field can be read in from the forcing data set.
     *
     * @param fieldName Name of the field to check
     * @return True if the field is valid, false otherwise
     */
    [[nodiscard]] bool isValid(const std::string& fieldName) const override
    {
        return forcingStateBefore.count(fieldName) && forcingStateAfter.count(fieldName);
    }

    // ModelComponent overrides
    std::string getName() const override { return "GeoreferencedInput"; };
    void setData(const ModelState::DataMap&) override { };

private:
    // Useful structs and aliases
    template <typename T> using RawDataMap = std::map<std::string, std::vector<T>>;
#ifdef USE_KOKKOS
    using ForcingDataMap = std::map<std::string, DeviceViewMA>;
#else
    using ForcingDataMap = ModelState::DataMap;
#endif

    struct {
        TimePoint before, after;
    } timeRange;

    // Weights and coordinate pointers for the bi-linear interpolation
    std::vector<double> xi, eta;
    std::vector<size_t> ij00, ij01, ij10, ij11;

    // Some (hopefully) self explanatory state variables
    TimePoint currentTime;
    ForcingDataMap forcingStateBefore, forcingStateAfter;
    std::string pathSpec, ncTimeName, ncLonName, ncLatName;
    ModelArray modelLons, modelLats;
    RawDataMap<double> forcingLonLats;
    std::set<std::string> forcings;
    std::set<std::pair<std::string, std::string>> vectors;
    std::unique_ptr<VectorRotator> rotator;
    std::vector<size_t> gridDims, gridStart;
    bool lonLat1D;

    // Read the netCDF file dimensions
    void readDims();

    // Reduce the read-in domain
    void tightenGrid();

    // Basic weight-setting functions
    void setWeights();
    void setWeights1D();
    void setWeights2D();

    /* Does a recursive search for the grid cell with corners i, ii, j, and jj, in which the point
     * {targetLon, targetLat} is found.
     */
    [[nodiscard]] bool recursiveBisectSearch(size_t k, FloatType targetLon, FloatType targetLat,
        size_t i, size_t ii, size_t j, size_t jj);

    /* Search for the grid cell in which the point {targetLon, targetLat} is found, by walking
     * through the source grid.
     */
    void findCellByWalking(FloatType targetLon, FloatType targetLat, size_t k);

    // Project {lon, lat} onto the orthographic coordinates {x,y} with {lon0, lat0} at its centre.
    static void orthographicProjection(
        double lon, double lat, FloatType lon0, FloatType lat0, double& x, double& y);

    /* Do a quick axis-aligned bounding box check to see if a point is (probably) in the grid cell.
     * This is a necessary condition, not a sufficient one.
     */
    [[nodiscard]] static bool pointInBoundingBox(
        const std::vector<double>& xCorners, const std::vector<double>& yCorners);

    /* Find the local coordinates and weights {xi, eta} for a bi-linear interpolation on a
     * curvilinear grid
     */
    [[nodiscard]] bool findLocalCoordinates(size_t k, double x00, double y00, double x10,
        double y10, double x01, double y01, double x11, double y11);

    // Apply the weights to do a bi-linear interpolation
    [[nodiscard]] ForcingDataMap interpolateSpatially(const RawDataMap<FloatType>& rawData) const;

    // Rotate the vectors from the input to model grid
    void rotateInputVectors(RawDataMap<FloatType>& rawData);

    // Read the forcing listed in ``forcings`` at times bracketing ``currentTime``.
    void readRawForcing(RawDataMap<FloatType>& rawDataBefore, RawDataMap<FloatType>& rawDataAfter);

    // The function that actually reads data from the netCDF file
    template <typename T>
    [[nodiscard]] RawDataMap<T> readRawData(
        const TimePoint& time, const std::set<std::string>& fields, size_t timeIndex = 0) const;

    // Wrap longitudes, so that lon0 is the largest value (usually either [-180 180] or [0 360]
    [[nodiscard]] static FloatType wrapLon(const double lon, const double lon0 = 0.)
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

    // Transform netCDF time to our time units
    [[nodiscard]] std::vector<TimePoint> getTimeFromNetCDF(const netCDF::NcFile& ncFile) const;
};

}

#endif // GEOREFERENCEDINPUT_HPP
