/*!
 * @file    Xios.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    31 July 2024
 * @brief   XIOS interface implementation
 * @details
 *
 * Implementation of XIOS interface
 *
 * This C++ interface is designed to implement core functionality of XIOS so
 * that it can be used in nextSIM-DG.
 *
 * To enable XIOS in nextSIM-DG add the following lines to the config file.
 *   [xios]
 *   enable = true
 */
#include <boost/date_time/posix_time/time_parsers.hpp>
#if USE_XIOS

#include "include/ModelArray.hpp"
#include "include/Xios.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <include/xios_c_interface.hpp>
#include <iostream>
#include <mpi.h>
#include <regex>
#include <string>

namespace Nextsim {

template <>
const std::map<int, std::string> Configured<Xios>::keyMap
    = { { Xios::ENABLED_KEY, "xios.enable" } };

/*!
 * Constructor
 *
 * Configure an XIOS server
 */
Xios::Xios() { configure(); }

//! Destructor
Xios::~Xios() { finalize(); }

//! Close XIOS context definition once xml config has been read and calendar settings updated
void Xios::close_context_definition()
{
    if (isEnabled) {
        cxios_context_close_definition();
    }
}

//! Finalize XIOS context
void Xios::context_finalize()
{
    if (isEnabled) {
        cxios_context_finalize();
    }
}

//! Close context and finialize server
void Xios::finalize()
{
    if (isEnabled) {
        cxios_finalize();
    }
}

/*!
 * Overrides `Configure` method from `Configured`
 *
 * Configure the XIOS server if XIOS is enabled in the settings.
 *
 */
void Xios::configure()
{
    // Check if XIOS is enabled in the nextSIM-DG configuration
    istringstream(Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()))
        >> std::boolalpha >> isEnabled;
    if (isEnabled) {
        configureServer();
    }
}

//! Configure calendar settings
void Xios::configureServer(const std::string calendarType)
{
    // Initialize XIOS Server process and store MPI communicator
    clientId = "client";
    nullComm_F = MPI_Comm_c2f(MPI_COMM_NULL);
    cxios_init_client(clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F);

    // Initialize MPI rank and size
    clientComm = MPI_Comm_f2c(clientComm_F);
    MPI_Comm_rank(clientComm, &mpi_rank);
    MPI_Comm_size(clientComm, &mpi_size);

    // Initialize 'nextSIM-DG' context
    contextId = "nextSIM-DG";
    cxios_context_initialize(contextId.c_str(), contextId.length(), &clientComm_F);

    // Initialize calendar wrapper for 'nextSIM-DG' context
    cxios_get_current_calendar_wrapper(&clientCalendar);
    cxios_set_calendar_wrapper_type(clientCalendar, calendarType.c_str(), calendarType.length());
    cxios_create_calendar(clientCalendar);
}

/*!
 * @return size of the client MPI communicator
 */
int Xios::getClientMPISize() { return mpi_size; }

/*!
 * @return rank of the client MPI communicator
 */
int Xios::getClientMPIRank() { return mpi_rank; }

/*!
 * Verify XIOS server is initialized
 *
 * @return true when XIOS server is initialized
 */
bool Xios::isInitialized()
{
    bool init { false };
    cxios_context_is_initialized(contextId.c_str(), contextId.length(), &init);
    return init;
}

/*!
 * Return datetime as std::string using ISO 8601 format (default).
 *
 * - If `isoFormat` is true  format will be 2023-03-03T17:11:00Z
 * - If `isoFormat` is false format will be 2023-03-03 17:11:00
 *
 * @param XIOS datetime representation
 * @param isoFormat as bool
 * @return corresponding string representation
 */
std::string Xios::convertXiosDatetimeToString(const cxios_date datetime, const bool isoFormat)
{
    boost::format fmt;
    if (isoFormat) {
        fmt = boost::format("%1$4d-%2$02d-%3$02dT%4$02d:%5$02d:%6$02dZ") % datetime.year
            % datetime.month % datetime.day % datetime.hour % datetime.minute % datetime.second;
    } else {
        fmt = boost::format("%1$4d-%2$02d-%3$02d %4$02d:%5$02d:%6$02d") % datetime.year
            % datetime.month % datetime.day % datetime.hour % datetime.minute % datetime.second;
    }
    return fmt.str();
}

/*!
 * Return std::string in ISO 8601 format (default) as an XIOS datetime object.
 *
 * - If `isoFormat` is true  format will be 2023-03-03T17:11:00Z
 * - If `isoFormat` is false format will be 2023-03-03 17:11:00
 *
 * @param string representation
 * @param isoFormat as bool
 * @return corresponding XIOS datetime representation
 */
cxios_date Xios::convertStringToXiosDatetime(const std::string datetimeStr, const bool isoFormat)
{
    std::string str = datetimeStr;
    if (isoFormat) {
        str = str.replace(10, 1, " "); // replaces T with a space
        str = str.replace(19, 1, " "); // replaces Z with a space
    }
    return cxios_date_convert_from_string(str.c_str(), str.length());
}

/*!
 * Set calendar origin
 *
 * @param origin
 */
void Xios::setCalendarOrigin(const TimePoint origin)
{
    cxios_date datetime = convertStringToXiosDatetime(origin.format(), true);
    cxios_set_calendar_wrapper_date_time_origin(clientCalendar, datetime);
}

/*!
 * Set calendar start date
 *
 * @param start date
 */
void Xios::setCalendarStart(const TimePoint start)
{
    cxios_date datetime = convertStringToXiosDatetime(start.format(), true);
    cxios_set_calendar_wrapper_date_start_date(clientCalendar, datetime);
}

/*!
 * Set calendar timestep
 *
 * @param timestep
 */
void Xios::setCalendarTimestep(const Duration timestep)
{
    cxios_duration duration { 0.0, 0.0, 0.0, 0.0, 0.0, timestep.seconds() };
    cxios_set_calendar_wrapper_timestep(clientCalendar, duration);
    cxios_update_calendar_timestep(clientCalendar);
}

/*!
 * Get calendar type
 *
 * @return calendar type
 */
std::string Xios::getCalendarType()
{
    char cStr[cStrLen];
    cxios_get_calendar_wrapper_type(clientCalendar, cStr, cStrLen);
    std::string calendarType(cStr, cStrLen);
    boost::algorithm::trim_right(calendarType);
    return calendarType;
}

/*!
 * Get calendar origin
 *
 * @return calendar origin
 */
TimePoint Xios::getCalendarOrigin()
{
    cxios_date calendar_origin;
    cxios_get_calendar_wrapper_date_time_origin(clientCalendar, &calendar_origin);
    return TimePoint(convertXiosDatetimeToString(calendar_origin, true));
}

/*!
 * Get calendar start date
 *
 * @return calendar start date
 */
TimePoint Xios::getCalendarStart()
{
    cxios_date calendar_start;
    cxios_get_calendar_wrapper_date_start_date(clientCalendar, &calendar_start);
    return TimePoint(convertXiosDatetimeToString(calendar_start, true));
}

/*!
 * Get calendar timestep
 *
 * @return calendar timestep
 */
Duration Xios::getCalendarTimestep()
{
    cxios_duration calendar_timestep;
    cxios_get_calendar_wrapper_timestep(clientCalendar, &calendar_timestep);
    char cStr[cStrLen];
    cxios_duration_convert_to_string(calendar_timestep, cStr, cStrLen);
    std::string durationStr(cStr, cStrLen);
    boost::algorithm::trim_right(durationStr);
    boost::erase_all(durationStr, "s");
    return Duration(std::stod(durationStr));
}

/*!
 * Get calendar step
 *
 * @return calendar step
 */
int Xios::getCalendarStep() { return clientCalendar->getCalendar()->getStep(); }

/*!
 * Get current calendar date
 *
 * @return current calendar date
 */
std::string Xios::getCurrentDate(const bool isoFormat)
{
    cxios_date xiosDate;
    cxios_get_current_date(&xiosDate);
    return convertXiosDatetimeToString(xiosDate, isoFormat);
}

/*!
 * Update XIOS calendar iteration/step number
 *
 * @param current step number
 */
void Xios::updateCalendar(const int stepNumber) { cxios_update_calendar(stepNumber); }

/*!
 * Get the axis_definition group
 *
 * @return a pointer to the XIOS CAxisGroup object
 */
xios::CAxisGroup* Xios::getAxisGroup()
{
    const std::string groupId = { "axis_definition" };
    xios::CAxisGroup* group = NULL;
    cxios_axisgroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for axis_definition group");
    }
    return group;
}

/*!
 * Get the axis associated with a given ID
 *
 * @param the axis ID
 * @return a pointer to the XIOS CAxis object
 */
xios::CAxis* Xios::getAxis(const std::string axisId)
{
    xios::CAxis* axis = NULL;
    cxios_axis_handle_create(&axis, axisId.c_str(), axisId.length());
    if (!axis) {
        throw std::runtime_error("Xios: Null pointer for axis with ID '" + axisId + "'");
    }
    return axis;
}

/*!
 * Create an axis with some ID
 *
 * @param the axis ID
 */
void Xios::createAxis(const std::string axisId)
{
    xios::CAxis* axis = NULL;
    cxios_xml_tree_add_axis(getAxisGroup(), &axis, axisId.c_str(), axisId.length());
    if (!axis) {
        throw std::runtime_error("Xios: Null pointer for axis with ID '" + axisId + "'");
    }
}

/*!
 * Set the size of a given axis (the number of global points)
 *
 * @param the axis ID
 * @param the size to set
 */
void Xios::setAxisSize(const std::string axisId, const size_t size)
{
    cxios_set_axis_n_glo(getAxis(axisId), (int)size);
}

/*!
 * Set the values associated with a given axis
 *
 * @param the axis ID
 * @param the values to set
 */
void Xios::setAxisValues(const std::string axisId, std::vector<double> values)
{
    if (!isDefinedAxisSize(axisId)) {
        setAxisSize(axisId, values.size());
    }
    int size = getAxisSize(axisId);
    if (size != values.size()) {
        throw std::runtime_error("Xios: Axis size incompatible with values for '" + axisId + "'");
    }
    cxios_set_axis_value(getAxis(axisId), values.data(), &size);
}

/*!
 * Get the size of a given axis (the number of global points)
 *
 * @param the axis ID
 * @return size of the corresponding axis
 */
size_t Xios::getAxisSize(const std::string axisId)
{
    int size;
    cxios_get_axis_n_glo(getAxis(axisId), &size);
    return (size_t)size;
}

/*!
 * Get the values associated with a given axis
 *
 * @param the axis ID
 * @return the corresponding values
 */
std::vector<double> Xios::getAxisValues(const std::string axisId)
{
    int size = getAxisSize(axisId);
    double* values = new double[size];
    cxios_get_axis_value(getAxis(axisId), values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
}

/*!
 * Verify whether a size has been defined for a given axis ID
 *
 * @param the axis ID
 * @return `true` if the size has been set, otherwise `false`
 */
bool Xios::isDefinedAxisSize(const std::string axisId)
{
    return cxios_is_defined_axis_n_glo(getAxis(axisId));
}

/*!
 * Verify whether values have been defined for a given axis ID
 *
 * @param the axis ID
 * @return `true` if the values have been set, otherwise `false`
 */
bool Xios::areDefinedAxisValues(const std::string axisId)
{
    return cxios_is_defined_axis_value(getAxis(axisId));
}

/*!
 * Get the domain_definition group
 *
 * @return a pointer to the XIOS CDomainGroup object
 */
xios::CDomainGroup* Xios::getDomainGroup()
{
    const std::string groupId = { "domain_definition" };
    xios::CDomainGroup* group = NULL;
    cxios_domaingroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for domain_definition group");
    }
    return group;
}

/*!
 * Get the domain associated with a given ID
 *
 * @param the domain ID
 * @return a pointer to the XIOS CDomain object
 */
xios::CDomain* Xios::getDomain(const std::string domainId)
{
    xios::CDomain* domain = NULL;
    cxios_domain_handle_create(&domain, domainId.c_str(), domainId.length());
    if (!domain) {
        throw std::runtime_error("Xios: Null pointer for domain with ID '" + domainId + "'");
    }
    return domain;
}

/*!
 * Create a domain with some ID
 *
 * @param the domain ID
 */
void Xios::createDomain(const std::string domainId)
{
    xios::CDomain* domain = NULL;
    cxios_xml_tree_add_domain(getDomainGroup(), &domain, domainId.c_str(), domainId.length());
    if (!domain) {
        throw std::runtime_error("Xios: Null pointer for domain with ID '" + domainId + "'");
    }
}

/*!
 * Set the local longitude size for a given domain
 *
 * @param the domain ID
 * @param the local longitude size
 */
void Xios::setDomainLongitudeSize(const std::string domainId, const size_t size)
{
    cxios_set_domain_ni(getDomain(domainId), (int)size);
}

/*!
 * Set the local latitude size for a given domain
 *
 * @param the domain ID
 * @param the local longitude size
 */
void Xios::setDomainLatitudeSize(const std::string domainId, const size_t size)
{
    cxios_set_domain_nj(getDomain(domainId), (int)size);
}

/*!
 * Set the local start longitude for a given domain
 *
 * @param the domain ID
 * @return the local start longitude
 */
void Xios::setDomainLongitudeStart(const std::string domainId, const size_t start)
{
    cxios_set_domain_ibegin(getDomain(domainId), (int)start);
}

/*!
 * Set the local start latitude for a given domain
 *
 * @param the domain ID
 * @return the local start latitude
 */
void Xios::setDomainLatitudeStart(const std::string domainId, const size_t start)
{
    cxios_set_domain_jbegin(getDomain(domainId), (int)start);
}

/*!
 * Set the local longitude values for a given domain
 *
 * @param the domain ID
 * @return the local longitude values
 */
void Xios::setDomainLongitudeValues(const std::string domainId, std::vector<double> values)
{
    int size = getDomainLongitudeSize(domainId);
    cxios_set_domain_lonvalue_1d(getDomain(domainId), values.data(), &size);
}

/*!
 * Set the local latitude values for a given domain
 *
 * @param the domain ID
 * @return the local latitude values
 */
void Xios::setDomainLatitudeValues(const std::string domainId, std::vector<double> values)
{
    int size = getDomainLatitudeSize(domainId);
    cxios_set_domain_latvalue_1d(getDomain(domainId), values.data(), &size);
}

/*!
 * Set the type of a given domain
 *
 * @param the domain ID
 * @param domain type to set
 */
void Xios::setDomainType(const std::string domainId, const std::string domainType)
{
    cxios_set_domain_type(getDomain(domainId), domainType.c_str(), domainType.length());
}

/*!
 * Set the global longitude size for a given domain
 *
 * @param the domain ID
 * @param global longitude size to set
 */
void Xios::setDomainGlobalLongitudeSize(const std::string domainId, const size_t size)
{
    cxios_set_domain_ni_glo(getDomain(domainId), (int)size);
}

/*!
 * Set the global latitude size for a given domain
 *
 * @param the domain ID
 * @param global latitude size to set
 */
void Xios::setDomainGlobalLatitudeSize(const std::string domainId, const size_t size)
{
    cxios_set_domain_nj_glo(getDomain(domainId), (int)size);
}

/*!
 * Get the type of a given domain
 *
 * @param the domain ID
 * @return the corresponding domain type
 */
std::string Xios::getDomainType(const std::string domainId)
{
    char cStr[cStrLen];
    cxios_get_domain_type(getDomain(domainId), cStr, cStrLen);
    std::string domainType(cStr, cStrLen);
    boost::algorithm::trim_right(domainType);
    return domainType;
}

/*!
 * Get the global longitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding global longitude size
 */
size_t Xios::getDomainGlobalLongitudeSize(const std::string domainId)
{
    int size;
    cxios_get_domain_ni_glo(getDomain(domainId), &size);
    return (size_t)size;
}

/*!
 * Get the global latitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding global latitude size
 */
size_t Xios::getDomainGlobalLatitudeSize(const std::string domainId)
{
    int size;
    cxios_get_domain_nj_glo(getDomain(domainId), &size);
    return (size_t)size;
}

/*!
 * Get the local longitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding local longitude size
 */
size_t Xios::getDomainLongitudeSize(const std::string domainId)
{
    int size;
    cxios_get_domain_ni(getDomain(domainId), &size);
    return (size_t)size;
}

/*!
 * Get the local latitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding local latitude size
 */
size_t Xios::getDomainLatitudeSize(const std::string domainId)
{
    int size;
    cxios_get_domain_nj(getDomain(domainId), &size);
    return (size_t)size;
}

/*!
 * Get the local starting longitude for a given domain
 *
 * @param the domain ID
 * @return the local starting longitude of the corresponding domain
 */
size_t Xios::getDomainLongitudeStart(const std::string domainId)
{
    int start;
    cxios_get_domain_ibegin(getDomain(domainId), &start);
    return (size_t)start;
}

/*!
 * Get the local starting latitude for a given domain
 *
 * @param the domain ID
 * @return the local starting latitude of the corresponding domain
 */
size_t Xios::getDomainLatitudeStart(const std::string domainId)
{
    int start;
    cxios_get_domain_jbegin(getDomain(domainId), &start);
    return (size_t)start;
}

/*!
 * Get the local longitude values for a given domain
 *
 * @param the domain ID
 * @return the local longitude values of the corresponding domain
 */
std::vector<double> Xios::getDomainLongitudeValues(const std::string domainId)
{
    int size = getDomainLongitudeSize(domainId);
    double* values = new double[size];
    cxios_get_domain_lonvalue_1d(getDomain(domainId), values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
}

/*!
 * Get the local latitude values for a given domain
 *
 * @param the domain ID
 * @return the local latitude values of the corresponding domain
 */
std::vector<double> Xios::getDomainLatitudeValues(const std::string domainId)
{
    int size = getDomainLatitudeSize(domainId);
    double* values = new double[size];
    cxios_get_domain_latvalue_1d(getDomain(domainId), values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
}

/*!
 * Verify whether a type has been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the type has been set, otherwise `false`
 */
bool Xios::isDefinedDomainType(const std::string domainId)
{
    return cxios_is_defined_domain_type(getDomain(domainId));
}

/*!
 * Verify whether a global longitude size has been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the global longitude size has been set, otherwise `false`
 */
bool Xios::isDefinedDomainGlobalLongitudeSize(const std::string domainId)
{
    return cxios_is_defined_domain_ni_glo(getDomain(domainId));
}

/*!
 * Verify whether a global latitude size has been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the global latitude size has been set, otherwise `false`
 */
bool Xios::isDefinedDomainGlobalLatitudeSize(const std::string domainId)
{
    return cxios_is_defined_domain_nj_glo(getDomain(domainId));
}

/*!
 * Verify whether a local longitude size has been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the local longitude size has been set, otherwise `false`
 */
bool Xios::isDefinedDomainLongitudeSize(const std::string domainId)
{
    return cxios_is_defined_domain_ni(getDomain(domainId));
}

/*!
 * Verify whether a local latitude size has been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the local latitude size has been set, otherwise `false`
 */
bool Xios::isDefinedDomainLatitudeSize(const std::string domainId)
{
    return cxios_is_defined_domain_nj(getDomain(domainId));
}

/*!
 * Verify whether a local starting longitude has been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the local starting longitude has been set, otherwise `false`
 */
bool Xios::isDefinedDomainLongitudeStart(const std::string domainId)
{
    return cxios_is_defined_domain_ibegin(getDomain(domainId));
}

/*!
 * Verify whether a local starting latitude has been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the local starting latitude has been set, otherwise `false`
 */
bool Xios::isDefinedDomainLatitudeStart(const std::string domainId)
{
    return cxios_is_defined_domain_jbegin(getDomain(domainId));
}

/*!
 * Verify whether a local longitude values have been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the local longitude values have been set, otherwise `false`
 */
bool Xios::areDefinedDomainLongitudeValues(const std::string domainId)
{
    return cxios_is_defined_domain_lonvalue_1d(getDomain(domainId));
}

/*!
 * Verify whether a local latitude values have been defined for a given domain ID
 *
 * @param the domain ID
 * @return `true` if the local latitude values have been set, otherwise `false`
 */
bool Xios::areDefinedDomainLatitudeValues(const std::string domainId)
{
    return cxios_is_defined_domain_latvalue_1d(getDomain(domainId));
}

/*!
 * Get the grid_definition group
 *
 * @return a pointer to the XIOS CGridGroup object
 */
xios::CGridGroup* Xios::getGridGroup()
{
    const std::string groupId = { "grid_definition" };
    xios::CGridGroup* group = NULL;
    cxios_gridgroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for grid_definition group");
    }
    return group;
}

/*!
 * Get the grid associated with a given ID
 *
 * @param the grid ID
 * @return a pointer to the XIOS CGrid object
 */
xios::CGrid* Xios::getGrid(const std::string gridId)
{
    xios::CGrid* grid = NULL;
    cxios_grid_handle_create(&grid, gridId.c_str(), gridId.length());
    if (!grid) {
        throw std::runtime_error("Xios: Null pointer for grid with ID '" + gridId + "'");
    }
    return grid;
}

/*!
 * Create a grid with some ID
 *
 * @param the grid ID
 */
void Xios::createGrid(const std::string gridId)
{
    xios::CGrid* grid = NULL;
    cxios_xml_tree_add_grid(getGridGroup(), &grid, gridId.c_str(), gridId.length());
    if (!grid) {
        throw std::runtime_error("Xios: Null pointer for grid with ID '" + gridId + "'");
    }
}

/*!
 * Get the name of a grid with a given ID
 *
 * @param the grid ID
 * @return name of the corresponding grid
 */
std::string Xios::getGridName(const std::string gridId)
{
    char cStr[cStrLen];
    cxios_get_grid_name(getGrid(gridId), cStr, cStrLen);
    std::string gridName(cStr, cStrLen);
    boost::algorithm::trim_right(gridName);
    return gridName;
}

/*!
 * Set the name of a grid with a given ID
 *
 * @param the grid ID
 * @param name to set
 */
void Xios::setGridName(const std::string gridId, const std::string gridName)
{
    cxios_set_grid_name(getGrid(gridId), gridName.c_str(), gridName.length());
}

/*!
 * Verify whether a name has been defined for a given grid ID
 *
 * @param the grid ID
 * @return `true` if the name has been set, otherwise `false`
 */
bool Xios::isDefinedGridName(const std::string gridId)
{
    return cxios_is_defined_grid_name(getGrid(gridId));
}

/*!
 * Associate an axis with a grid
 *
 * @param the grid ID
 * @param the axis ID
 */
void Xios::gridAddAxis(const std::string gridId, const std::string axisId)
{
    xios::CAxis* axis = getAxis(axisId);
    cxios_xml_tree_add_axistogrid(getGrid(gridId), &axis, axisId.c_str(), axisId.length());
}

/*!
 * Associate a domain with a grid
 *
 * @param the grid ID
 * @param the domain ID
 */
void Xios::gridAddDomain(const std::string gridId, const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    cxios_xml_tree_add_domaintogrid(getGrid(gridId), &domain, domainId.c_str(), domainId.length());
}

/*!
 * Get all axis IDs associated with a given grid
 *
 * @param the grid ID
 * @return all axis IDs associated with the grid
 */
std::vector<std::string> Xios::gridGetAxisIds(const std::string gridId)
{
    return getGrid(gridId)->getAxisList();
}

/*!
 * Get all domain IDs associated with a given grid
 *
 * @param the grid ID
 * @return all domain IDs associated with the grid
 */
std::vector<std::string> Xios::gridGetDomainIds(const std::string gridId)
{
    return getGrid(gridId)->getDomainList();
}

/*!
 * Get the field_definition group
 *
 * @return a pointer to the XIOS CFieldGroup object
 */
xios::CFieldGroup* Xios::getFieldGroup()
{
    const std::string groupId = { "field_definition" };
    xios::CFieldGroup* group = NULL;
    cxios_fieldgroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for field_definition group");
    }
    return group;
}

/*!
 * Get the field associated with a given ID
 *
 * @param the field ID
 * @return a pointer to the XIOS CField object
 */
xios::CField* Xios::getField(const std::string fieldId)
{
    xios::CField* field = NULL;
    cxios_field_handle_create(&field, fieldId.c_str(), fieldId.length());
    if (!field) {
        throw std::runtime_error("Xios: Null pointer for field with ID '" + fieldId + "'");
    }
    return field;
}

/*!
 * Create a field with some ID
 *
 * @param the field ID
 */
void Xios::createField(const std::string fieldId)
{
    xios::CField* field = NULL;
    cxios_xml_tree_add_field(getFieldGroup(), &field, fieldId.c_str(), fieldId.length());
    if (!field) {
        throw std::runtime_error("Xios: Null pointer for field with ID '" + fieldId + "'");
    }
}

/*!
 * Set the name of a field with a given ID
 *
 * @param the field ID
 * @param name to set
 */
void Xios::setFieldName(const std::string fieldId, const std::string fieldName)
{
    cxios_set_field_name(getField(fieldId), fieldName.c_str(), fieldName.length());
}

/*!
 * Set the operation for a field with a given ID
 *
 * @param the field ID
 * @param operation to set
 */
void Xios::setFieldOperation(const std::string fieldId, const std::string operation)
{
    cxios_set_field_operation(getField(fieldId), operation.c_str(), operation.length());
}

/*!
 * Set the grid reference for a field with a given ID
 *
 * @param the field ID
 * @param grid reference to set
 */
void Xios::setFieldGridRef(const std::string fieldId, const std::string gridRef)
{
    cxios_set_field_grid_ref(getField(fieldId), gridRef.c_str(), gridRef.length());
}

/*!
 * Get the name of a field with a given ID
 *
 * @param the field ID
 * @return name of the corresponding field
 */
std::string Xios::getFieldName(const std::string fieldId)
{
    char cStr[cStrLen];
    cxios_get_field_name(getField(fieldId), cStr, cStrLen);
    std::string fieldName(cStr, cStrLen);
    boost::algorithm::trim_right(fieldName);
    return fieldName;
}

/*!
 * Get the operation associated with a field with a given ID
 *
 * @param the field ID
 * @return operation used for the corresponding field
 */
std::string Xios::getFieldOperation(const std::string fieldId)
{
    char cStr[cStrLen];
    cxios_get_field_operation(getField(fieldId), cStr, cStrLen);
    std::string operation(cStr, cStrLen);
    boost::algorithm::trim_right(operation);
    return operation;
}

/*!
 * Get the grid reference associated with a field with a given ID
 *
 * @param the field ID
 * @return grid reference used for the corresponding field
 */
std::string Xios::getFieldGridRef(const std::string fieldId)
{
    char cStr[cStrLen];
    cxios_get_field_grid_ref(getField(fieldId), cStr, cStrLen);
    std::string gridRef(cStr, cStrLen);
    boost::algorithm::trim_right(gridRef);
    return gridRef;
}

/*!
 * Verify whether a name has been defined for a given field ID
 *
 * @param the field ID
 * @return `true` if the name has been set, otherwise `false`
 */
bool Xios::isDefinedFieldName(const std::string fieldId)
{
    return cxios_is_defined_field_name(getField(fieldId));
}

/*!
 * Verify whether an operation has been defined for a given field ID
 *
 * @param the field ID
 * @return `true` if the operation has been set, otherwise `false`
 */
bool Xios::isDefinedFieldOperation(const std::string fieldId)
{
    return cxios_is_defined_field_operation(getField(fieldId));
}

/*!
 * Verify whether a grid reference has been defined for a given field ID
 *
 * @param the field ID
 * @return `true` if the grid reference has been set, otherwise `false`
 */
bool Xios::isDefinedFieldGridRef(const std::string fieldId)
{
    return cxios_is_defined_field_grid_ref(getField(fieldId));
}

/*!
 * Get the file_definition group
 *
 * @return a pointer to the XIOS CFileGroup object
 */
xios::CFileGroup* Xios::getFileGroup()
{
    const std::string groupId = { "file_definition" };
    xios::CFileGroup* group = NULL;
    cxios_filegroup_handle_create(&group, groupId.c_str(), groupId.length());
    if (!group) {
        throw std::runtime_error("Xios: Null pointer for file_definition group");
    }
    return group;
}

/*!
 * Get the file associated with a given ID
 *
 * @param the file ID
 * @return a pointer to the XIOS CFile object
 */
xios::CFile* Xios::getFile(const std::string fileId)
{
    xios::CFile* file = NULL;
    cxios_file_handle_create(&file, fileId.c_str(), fileId.length());
    if (!file) {
        throw std::runtime_error("Xios: Null pointer for file with ID '" + fileId + "'");
    }
    return file;
}

/*!
 * Create a file with some ID
 *
 * @param the file ID
 */
void Xios::createFile(const std::string fileId)
{
    xios::CFile* file = NULL;
    cxios_xml_tree_add_file(getFileGroup(), &file, fileId.c_str(), fileId.length());
    if (!file) {
        throw std::runtime_error("Xios: Null pointer for file with ID '" + fileId + "'");
    }
}

/*!
 * Set the name of a file with a given ID
 *
 * @param the file ID
 * @param file name to set
 */
void Xios::setFileName(const std::string fileId, const std::string fileName)
{
    cxios_set_file_name(getFile(fileId), fileName.c_str(), fileName.length());
}

/*!
 * Set the type of a file with a given ID
 *
 * @param the file ID
 * @param file type to set
 */
void Xios::setFileType(const std::string fileId, const std::string fileType)
{
    cxios_set_file_type(getFile(fileId), fileType.c_str(), fileType.length());
}

/*!
 * Set the output frequency of a file with a given ID
 *
 * @param the file ID
 * @param output frequency to set
 */
void Xios::setFileOutputFreq(const std::string fileId, const std::string freq)
{
    cxios_set_file_output_freq(
        getFile(fileId), cxios_duration_convert_from_string(freq.c_str(), freq.length()));
}

/*!
 * Set the split frequency of a file with a given ID
 *
 * @param the file ID
 * @param split frequency to set
 */
void Xios::setFileSplitFreq(const std::string fileId, const std::string freq)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_split_freq(file)) {
        Logged::warning("Xios: Split frequency already set for file '" + fileId + "'");
    }
    cxios_set_file_split_freq(
        file, cxios_duration_convert_from_string(freq.c_str(), freq.length()));
    if (!cxios_is_defined_file_split_freq(file)) {
        throw std::runtime_error("Xios: Failed to set split frequency for file '" + fileId + "'");
    }
}

/*!
 * Get the name of a file with a given ID
 *
 * @param the file ID
 * @return name of the corresponding file
 */
std::string Xios::getFileName(const std::string fileId)
{
    char cStr[cStrLen];
    cxios_get_file_name(getFile(fileId), cStr, cStrLen);
    std::string fileName(cStr, cStrLen);
    boost::algorithm::trim_right(fileName);
    return fileName;
}

/*!
 * Get the type of a file with a given ID
 *
 * @param the file ID
 * @return type of the corresponding file
 */
std::string Xios::getFileType(const std::string fileId)
{
    char cStr[cStrLen];
    cxios_get_file_type(getFile(fileId), cStr, cStrLen);
    std::string fileType(cStr, cStrLen);
    boost::algorithm::trim_right(fileType);
    return fileType;
}

/*!
 * Get the output frequency of a file with a given ID
 *
 * @param the file ID
 * @return the corresponding output frequency
 */
std::string Xios::getFileOutputFreq(const std::string fileId)
{
    cxios_duration duration;
    cxios_get_file_output_freq(getFile(fileId), &duration);
    char cStr[cStrLen];
    cxios_duration_convert_to_string(duration, cStr, cStrLen);
    std::string outputFreq(cStr, cStrLen);
    boost::algorithm::trim_right(outputFreq);
    return outputFreq;
}

/*!
 * Get the split frequency of a file with a given ID
 *
 * @param the file ID
 * @return split frequency of the corresponding file
 */
std::string Xios::getFileSplitFreq(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_split_freq(file)) {
        throw std::runtime_error("Xios: Undefined values for file '" + fileId + "'");
    }
    cxios_duration duration;
    cxios_get_file_split_freq(file, &duration);
    char cStr[cStrLen];
    cxios_duration_convert_to_string(duration, cStr, cStrLen);
    std::string freq(cStr, cStrLen);
    boost::algorithm::trim_right(freq);
    return freq;
}

/*!
 * Verify whether a given file ID is valid
 *
 * @param the file ID
 * @return `true` if the file ID is valid, otherwise `false`
 */
bool Xios::validFileId(const std::string fileId)
{
    bool valid;
    cxios_file_valid_id(&valid, fileId.c_str(), fileId.length());
    return valid;
}

/*!
 * Get all field IDs associated with a given file
 *
 * @param the file ID
 * @return all field IDs associated with the file
 */
std::vector<std::string> Xios::fileGetFieldIds(const std::string fileId)
{
    std::vector<xios::CField*> fields = getFile(fileId)->getAllFields();
    std::vector<std::string> fieldIds(fields.size());
    for (int i = 0; i < fields.size(); i++) {
        fieldIds[i] = fields[i]->getId();
    }
    return fieldIds;
}

/*!
 * Verify whether a name has been defined for a given file ID
 *
 * @param the file ID
 * @return `true` if the name has been set, otherwise `false`
 */
bool Xios::isDefinedFileName(const std::string fileId)
{
    return cxios_is_defined_file_name(getFile(fileId));
}

/*!
 * Verify whether a type has been defined for a given file ID
 *
 * @param the file ID
 * @return `true` if the type has been set, otherwise `false`
 */
bool Xios::isDefinedFileType(const std::string fileId)
{
    return cxios_is_defined_file_type(getFile(fileId));
}

/*!
 * Verify whether an output frequency has been defined for a given file ID
 *
 * @param the file ID
 * @return `true` if the output frequency has been set, otherwise `false`
 */
bool Xios::isDefinedFileOutputFreq(const std::string fileId)
{
    return cxios_is_defined_file_output_freq(getFile(fileId));
}

/*!
 * Associate a field with a file
 *
 * @param the file ID
 * @param the field ID
 */
void Xios::fileAddField(const std::string fileId, const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    cxios_xml_tree_add_fieldtofile(getFile(fileId), &field, fieldId.c_str(), fieldId.length());
}

/*!
 * Send a field to the XIOS server to be written to file
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void Xios::write(const std::string fieldId, ModelArray& modelarray)
{
    auto ndim = modelarray.nDimensions();
    auto dims = modelarray.dimensions();
    if (ndim == 2) {
        cxios_write_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], -1);
    } else if (ndim == 3) {
        cxios_write_data_k83(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], dims[2], -1);
    } else if (ndim == 4) {
        cxios_write_data_k84(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], dims[2], dims[3], -1);
    } else {
        throw std::invalid_argument("Only ModelArrays of dimension 2, 3, or 4 are supported");
    }
}
}

#endif
