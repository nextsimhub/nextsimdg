/*!
 * @file NetCDFForcings.hpp
 *
 * @date Oct 21, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelArray.hpp"
#include "include/Time.hpp"

#include <functional>
#include <string>
#include <tuple>

#ifndef NETCDFFORCINGS_HPP
#define NETCDFFORCINGS_HPP

namespace Nextsim {

class NetCDFForcings {
public:
    using CoordinateIndexFn = std::function<std::pair<double, double>(double, double)>;
    using NameLookupFn = std::function<std::string(const std::string&)>;
    using FileLookupFn = std::function<std::string(const std::string&, const TimePoint&)>;
    using Buffer = Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

    /*!
     * @brief Sets the function which performs the look-up from coordinate
     * values to (fractional) indices on the source grid.
     *
     * @details The function takes longitude and latitude as two arguments.
     * The order and units do not matter as long as the supplied coordinates
     * match up with the data provided to the targetLon and targetLat arrays.
     *
     * @param fn A std::function that translates from two double values
     * longitude and latitude) to the x and y indices of the corresponding
     * point on the source data grid.
     */
    void setCoordinateIndexFn(const CoordinateIndexFn& fn);
    /*!
     * @brief Sets the function that maps nextSIM field names to field name in
     * the forcing file.
     *
     * @param fn A std::function that maps a string nextSIM field name to the
     *     name of the field in the forcing file.
     */
    void setFieldNameLookupFn(const NameLookupFn& fn);

    /*!
     * @brief Sets the function that maps from nextSIM field name and
     * TimePoint of interest to the filename that will contain the data.
     *
     * @param fn A std::function that maps a std::string and a TimePoint to a
     *     filename that contains that the corresponding data.
     */
    void setFileNameLookupFn(const FileLookupFn& fn);
    /*!
     * @brief Sets the longitude values of the grid onto which the forcing
     * data should be interpolated.
     *
     * @details The units (degrees or radians) and even whether this array
     * represents longitude at all are only required to match the inputs of the
     * function passed to NetCDFForcings::setCoordinateIndexFn(). The data in
     * the array set by this function will be passed as to coordinate index
     * calculation function as the first argument.
     *
     * @param longitudeArray. The first coordinate for the interpolation target
     *     grid.
     */
    void setTargetLongitude(const ModelArray& longitudeArray);

    /*!
     * @brief Sets the latitude values of the grid onto which the forcing
     * data should be interpolated.
     *
     * @details The units (degrees or radians) and even whether this array
     * represents latitude at all are only required to match the inputs of the
     * function passed to NetCDFForcings::setCoordinateIndexFn(). The data in
     * the array set by this function will be passed as to coordinate index
     * calculation function as the second argument.
     *
     * @param latitudeArray. The second coordinate for the interpolation target
     *     grid.
     */
    void setTargetLatitude(const ModelArray& latitudeArray);

    ModelArray getData(const std::string& nsName, const TimePoint& time);
    static Buffer getFileIndexData(const std::string& filename, const std::string& fieldName, size_t tIndex);
    static Buffer getFileIndexData(const std::string& filename, size_t tIndex);
    static Buffer getFileIndexData(const std::string& filename, const std::string& fieldName, size_t tIndex, double& missing);
    static Buffer getFileIndexData(const std::string& filename, size_t tIndex, double& missing);

    static ModelArray maFromBuffer(const Buffer& buffer, const ModelArray& fracI, const ModelArray& fracJ);
    static ModelArray maFromBuffer(const Buffer& buffer, const ModelArray& fracI, const ModelArray& fracJ, double missing);
private:

    const Buffer getBufferData(const std::string nsName, const TimePoint& time);

    CoordinateIndexFn indicesFromLonLat;
    NameLookupFn forcingNameFromNSname;
    FileLookupFn fileNameFn;

    ModelArray targetLon;
    ModelArray targetLat;
};

} /* namespace Nextsim */

#endif /* NETCDFFORCINGS_HPP */
