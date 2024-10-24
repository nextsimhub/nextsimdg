/*!
 * @file    Xios.cpp
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    21 August 2024
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

static const std::map<int, std::string> keyMap = { { Xios::ENABLED_KEY, "xios.enable" } };

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
 * Convert a C-string to a C++ `std::string`.
 *
 * @param C-string
 * @param length of C-string
 * @return C++ string version
 */
std::string Xios::convertCStrToCppStr(const char* cStr, int cStrLen)
{
    std::string cppStr(cStr, cStrLen);
    boost::algorithm::trim_right(cppStr);
    return cppStr;
}

/*!
 * Convert an XIOS duration object into a nextSIM-DG one.
 *
 * @param XIOS duration object
 * @return nextSIM-DG version
 */
Duration Xios::convertDurationFromXios(const cxios_duration duration)
{
    char cStr[cStrLen];
    cxios_duration_convert_to_string(duration, cStr, cStrLen);
    std::string durationStr = convertCStrToCppStr(cStr, cStrLen);
    boost::erase_all(durationStr, "s");
    return Duration(std::stod(durationStr));
}

/*!
 * Convert a nextSIM-DG duration object into an XIOS one.
 *
 * @param nextSIM-DG duration object
 * @return XIOS version
 */
cxios_duration Xios::convertDurationToXios(const Duration duration)
{
    return cxios_duration({ 0.0, 0.0, 0.0, 0.0, 0.0, duration.seconds() });
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
    cxios_set_calendar_wrapper_timestep(clientCalendar, convertDurationToXios(timestep));
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
    return convertCStrToCppStr(cStr, cStrLen);
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
    return convertDurationFromXios(calendar_timestep);
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
        throw std::runtime_error("Xios: Null pointer for group 'axis_definition'");
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
    bool exists;
    cxios_axis_valid_id(&exists, axisId.c_str(), axisId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined axis '" + axisId + "'");
    }
    xios::CAxis* axis = NULL;
    cxios_axis_handle_create(&axis, axisId.c_str(), axisId.length());
    if (!axis) {
        throw std::runtime_error("Xios: Null pointer for axis '" + axisId + "'");
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
    bool exists;
    cxios_axis_valid_id(&exists, axisId.c_str(), axisId.length());
    if (exists) {
        throw std::runtime_error("Xios: Axis '" + axisId + "' already exists");
    }
    xios::CAxis* axis = NULL;
    cxios_xml_tree_add_axis(getAxisGroup(), &axis, axisId.c_str(), axisId.length());
    if (!axis) {
        throw std::runtime_error("Xios: Null pointer for axis '" + axisId + "'");
    }
    cxios_axis_valid_id(&exists, axisId.c_str(), axisId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create axis '" + axisId + "'");
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
    xios::CAxis* axis = getAxis(axisId);
    if (cxios_is_defined_axis_n_glo(axis)) {
        Logged::warning("Xios: Size already set for axis '" + axisId + "'");
    }
    cxios_set_axis_n_glo(axis, (int)size);
    if (!cxios_is_defined_axis_n_glo(axis)) {
        throw std::runtime_error("Xios: Failed to set size for axis '" + axisId + "'");
    }
}

/*!
 * Set the values associated with a given axis
 *
 * @param the axis ID
 * @param the values to set
 */
void Xios::setAxisValues(const std::string axisId, std::vector<double> values)
{
    xios::CAxis* axis = getAxis(axisId);
    if (cxios_is_defined_axis_value(axis)) {
        Logged::warning("Xios: Values already set for axis '" + axisId + "'");
    }
    if (!cxios_is_defined_axis_n_glo(axis)) {
        setAxisSize(axisId, values.size());
    }
    int size = getAxisSize(axisId);
    if (size != values.size()) {
        throw std::runtime_error("Xios: Size incompatible with values for axis '" + axisId + "'");
    }
    cxios_set_axis_value(axis, values.data(), &size);
    if (!cxios_is_defined_axis_value(axis)) {
        throw std::runtime_error("Xios: Failed to set values for axis '" + axisId + "'");
    }
}

/*!
 * Get the size of a given axis (the number of global points)
 *
 * @param the axis ID
 * @return size of the corresponding axis
 */
size_t Xios::getAxisSize(const std::string axisId)
{
    xios::CAxis* axis = getAxis(axisId);
    if (!cxios_is_defined_axis_n_glo(axis)) {
        throw std::runtime_error("Xios: Undefined size for axis '" + axisId + "'");
    }
    int size;
    cxios_get_axis_n_glo(axis, &size);
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
    xios::CAxis* axis = getAxis(axisId);
    if (!cxios_is_defined_axis_value(axis)) {
        throw std::runtime_error("Xios: Undefined values for axis '" + axisId + "'");
    }
    int size = getAxisSize(axisId);
    double* values = new double[size];
    cxios_get_axis_value(axis, values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
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
        throw std::runtime_error("Xios: Null pointer for group 'domain_definition'");
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
    bool exists;
    cxios_domain_valid_id(&exists, domainId.c_str(), domainId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined domain '" + domainId + "'");
    }
    xios::CDomain* domain = NULL;
    cxios_domain_handle_create(&domain, domainId.c_str(), domainId.length());
    if (!domain) {
        throw std::runtime_error("Xios: Null pointer for domain '" + domainId + "'");
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
    bool exists;
    cxios_domain_valid_id(&exists, domainId.c_str(), domainId.length());
    if (exists) {
        throw std::runtime_error("Xios: Domain '" + domainId + "' already exists");
    }
    xios::CDomain* domain = NULL;
    cxios_xml_tree_add_domain(getDomainGroup(), &domain, domainId.c_str(), domainId.length());
    if (!domain) {
        throw std::runtime_error("Xios: Null pointer for domain '" + domainId + "'");
    }
    cxios_domain_valid_id(&exists, domainId.c_str(), domainId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create domain '" + domainId + "'");
    }
}

/*!
 * Set the global number of points in the x-direction for a given domain
 *
 * @param the domain ID
 * @param the global number of points in the x-direction
 */
void Xios::setDomainGlobalXSize(const std::string domainId, const size_t size)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_ni_glo(domain)) {
        Logged::warning("Xios: Overwriting global x-size for domain '" + domainId + "'");
    }
    cxios_set_domain_ni_glo(domain, (int)size);
    if (!cxios_is_defined_domain_ni_glo(domain)) {
        throw std::runtime_error("Xios: Failed to set global x-size for domain '" + domainId + "'");
    }
}

/*!
 * Set the global number of points in the y-direction for a given domain
 *
 * @param the domain ID
 * @param the global number of points in the y-direction
 */
void Xios::setDomainGlobalYSize(const std::string domainId, const size_t size)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_nj_glo(domain)) {
        Logged::warning("Xios: Overwriting global y-size for domain '" + domainId + "'");
    }
    cxios_set_domain_nj_glo(domain, (int)size);
    if (!cxios_is_defined_domain_nj_glo(domain)) {
        throw std::runtime_error("Xios: Failed to set global y-size for domain '" + domainId + "'");
    }
}

/*!
 * Set the local number of points in the x-direction for a given domain
 *
 * @param the domain ID
 * @param the local number of points in the x-direction
 */
void Xios::setDomainLocalXSize(const std::string domainId, const size_t size)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_ni(domain)) {
        Logged::warning("Xios: Overwriting local x-size for domain '" + domainId + "'");
    }
    cxios_set_domain_ni(domain, (int)size);
    if (!cxios_is_defined_domain_ni(domain)) {
        throw std::runtime_error("Xios: Failed to set local x-size for domain '" + domainId + "'");
    }
}

/*!
 * Set the local number of points in the y-direction for a given domain
 *
 * @param the domain ID
 * @param the local number of points in the y-direction
 */
void Xios::setDomainLocalYSize(const std::string domainId, const size_t size)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_nj(domain)) {
        Logged::warning("Xios: Overwriting local y-size for domain '" + domainId + "'");
    }
    cxios_set_domain_nj(domain, (int)size);
    if (!cxios_is_defined_domain_nj(domain)) {
        throw std::runtime_error("Xios: Failed to set local y-size for domain '" + domainId + "'");
    }
}

/*!
 * Set the local starting x-index for a given domain
 *
 * @param the domain ID
 * @return the local starting x-index
 */
void Xios::setDomainLocalXStart(const std::string domainId, const size_t start)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_ibegin(domain)) {
        Logged::warning("Xios: Overwriting local starting x-index for domain '" + domainId + "'");
    }
    cxios_set_domain_ibegin(domain, (int)start);
    if (!cxios_is_defined_domain_ibegin(domain)) {
        throw std::runtime_error(
            "Xios: Failed to set local starting x-index for domain '" + domainId + "'");
    }
}

/*!
 * Set the local starting y-index for a given domain
 *
 * @param the domain ID
 * @return the local starting y-index
 */
void Xios::setDomainLocalYStart(const std::string domainId, const size_t start)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_jbegin(domain)) {
        Logged::warning("Xios: Overwriting local starting y-index for domain '" + domainId + "'");
    }
    cxios_set_domain_jbegin(domain, (int)start);
    if (!cxios_is_defined_domain_jbegin(domain)) {
        throw std::runtime_error(
            "Xios: Failed to set local starting y-index for domain '" + domainId + "'");
    }
}

/*!
 * Set the local x-values for a given domain
 *
 * @param the domain ID
 * @return the local x-values
 */
void Xios::setDomainLocalXValues(const std::string domainId, std::vector<double> values)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_lonvalue_1d(domain)) {
        Logged::warning("Xios: Overwriting local x-values for domain '" + domainId + "'");
    }
    if (!cxios_is_defined_domain_ni(domain)) {
        setDomainLocalXSize(domainId, values.size());
    }
    int size = getDomainLocalXSize(domainId);
    cxios_set_domain_lonvalue_1d(domain, values.data(), &size);
    if (!cxios_is_defined_domain_lonvalue_1d(domain)) {
        throw std::runtime_error(
            "Xios: Failed to set local x-values for domain '" + domainId + "'");
    }
}

/*!
 * Set the local y-values for a given domain
 *
 * @param the domain ID
 * @return the local y-values
 */
void Xios::setDomainLocalYValues(const std::string domainId, std::vector<double> values)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_latvalue_1d(domain)) {
        Logged::warning("Xios: Overwriting local y-values for domain '" + domainId + "'");
    }
    if (!cxios_is_defined_domain_nj(domain)) {
        setDomainLocalYSize(domainId, values.size());
    }
    int size = getDomainLocalYSize(domainId);
    cxios_set_domain_latvalue_1d(domain, values.data(), &size);
    if (!cxios_is_defined_domain_latvalue_1d(domain)) {
        throw std::runtime_error(
            "Xios: Failed to set local y-values for domain '" + domainId + "'");
    }
}

/*!
 * Set the type of a given domain
 *
 * @param the domain ID
 * @param domain type to set
 */
void Xios::setDomainType(const std::string domainId, const std::string domainType)
{
    xios::CDomain* domain = getDomain(domainId);
    if (cxios_is_defined_domain_type(domain)) {
        Logged::warning("Xios: Overwriting type for domain '" + domainId + "'");
    }
    cxios_set_domain_type(domain, domainType.c_str(), domainType.length());
    if (!cxios_is_defined_domain_type(domain)) {
        throw std::runtime_error("Xios: Failed to set type for domain '" + domainId + "'");
    }
}

/*!
 * Get the type of a given domain
 *
 * @param the domain ID
 * @return the corresponding domain type
 */
std::string Xios::getDomainType(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_type(domain)) {
        throw std::runtime_error("Xios: Undefined type for domain '" + domainId + "'");
    }
    char cStr[cStrLen];
    cxios_get_domain_type(domain, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the global number of points in the x-direction for a given domain
 *
 * @param the domain ID
 * @return the global number of points in the x-direction
 */
size_t Xios::getDomainGlobalXSize(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_ni_glo(domain)) {
        throw std::runtime_error("Xios: Undefined global x-size for domain '" + domainId + "'");
    }
    int size;
    cxios_get_domain_ni_glo(domain, &size);
    return (size_t)size;
}

/*!
 * Get the global number of points in the y-direction for a given domain
 *
 * @param the domain ID
 * @return the global number of points in the y-direction
 */
size_t Xios::getDomainGlobalYSize(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_nj_glo(domain)) {
        throw std::runtime_error("Xios: Undefined global y-size for domain '" + domainId + "'");
    }
    int size;
    cxios_get_domain_nj_glo(domain, &size);
    return (size_t)size;
}

/*!
 * Get the local number of points in the x-direction for a given domain
 *
 * @param the domain ID
 * @return the local number of points in the x-direction
 */
size_t Xios::getDomainLocalXSize(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_ni(domain)) {
        throw std::runtime_error("Xios: Undefined local x-size for domain '" + domainId + "'");
    }
    int size;
    cxios_get_domain_ni(domain, &size);
    return (size_t)size;
}

/*!
 * Get the local number of points in the y-direction for a given domain
 *
 * @param the domain ID
 * @return the local number of points in the y-direction
 */
size_t Xios::getDomainLocalYSize(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_nj(domain)) {
        throw std::runtime_error("Xios: Undefined local y-size for domain '" + domainId + "'");
    }
    int size;
    cxios_get_domain_nj(domain, &size);
    return (size_t)size;
}

/*!
 * Get the local starting x-index for a given domain
 *
 * @param the domain ID
 * @return the local starting x-index
 */
size_t Xios::getDomainLocalXStart(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_ibegin(domain)) {
        throw std::runtime_error(
            "Xios: Undefined local starting x-index for domain '" + domainId + "'");
    }
    int start;
    cxios_get_domain_ibegin(domain, &start);
    return (size_t)start;
}

/*!
 * Get the local starting y-index for a given domain
 *
 * @param the domain ID
 * @return the local starting y-index
 */
size_t Xios::getDomainLocalYStart(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_jbegin(domain)) {
        throw std::runtime_error(
            "Xios: Undefined local starting y-index for domain '" + domainId + "'");
    }
    int start;
    cxios_get_domain_jbegin(domain, &start);
    return (size_t)start;
}

/*!
 * Get the local x-values for a given domain
 *
 * @param the domain ID
 * @return the local x-values
 */
std::vector<double> Xios::getDomainLocalXValues(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_lonvalue_1d(domain)) {
        throw std::runtime_error("Xios: Undefined local x-values for domain '" + domainId + "'");
    }
    int size = getDomainLocalXSize(domainId);
    double* values = new double[size];
    cxios_get_domain_lonvalue_1d(domain, values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
}

/*!
 * Get the local y-values for a given domain
 *
 * @param the domain ID
 * @return the local y-values
 */
std::vector<double> Xios::getDomainLocalYValues(const std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    if (!cxios_is_defined_domain_latvalue_1d(domain)) {
        throw std::runtime_error("Xios: Undefined local y-values for domain '" + domainId + "'");
    }
    int size = getDomainLocalYSize(domainId);
    double* values = new double[size];
    cxios_get_domain_latvalue_1d(domain, values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
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
        throw std::runtime_error("Xios: Null pointer for group 'grid_definition'");
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
    bool exists;
    cxios_grid_valid_id(&exists, gridId.c_str(), gridId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined grid '" + gridId + "'");
    }
    xios::CGrid* grid = NULL;
    cxios_grid_handle_create(&grid, gridId.c_str(), gridId.length());
    if (!grid) {
        throw std::runtime_error("Xios: Null pointer for grid '" + gridId + "'");
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
    bool exists;
    cxios_grid_valid_id(&exists, gridId.c_str(), gridId.length());
    if (exists) {
        throw std::runtime_error("Xios: Grid '" + gridId + "' already exists");
    }
    xios::CGrid* grid = NULL;
    cxios_xml_tree_add_grid(getGridGroup(), &grid, gridId.c_str(), gridId.length());
    if (!grid) {
        throw std::runtime_error("Xios: Null pointer for grid '" + gridId + "'");
    }
    cxios_grid_valid_id(&exists, gridId.c_str(), gridId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create grid '" + gridId + "'");
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
    xios::CGrid* grid = getGrid(gridId);
    if (!cxios_is_defined_grid_name(grid)) {
        throw std::runtime_error("Xios: Undefined name for grid '" + gridId + "'");
    }
    char cStr[cStrLen];
    cxios_get_grid_name(grid, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Set the name of a grid with a given ID
 *
 * @param the grid ID
 * @param name to set
 */
void Xios::setGridName(const std::string gridId, const std::string gridName)
{
    xios::CGrid* grid = getGrid(gridId);
    if (cxios_is_defined_grid_name(grid)) {
        Logged::warning("Xios: Overwriting name for grid '" + gridId + "'");
    }
    cxios_set_grid_name(grid, gridName.c_str(), gridName.length());
    if (!cxios_is_defined_grid_name(grid)) {
        throw std::runtime_error("Xios: Failed to set name for grid '" + gridId + "'");
    }
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
        throw std::runtime_error("Xios: Null pointer for group 'field_definition'");
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
    bool exists;
    cxios_field_valid_id(&exists, fieldId.c_str(), fieldId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined field '" + fieldId + "'");
    }
    xios::CField* field = NULL;
    cxios_field_handle_create(&field, fieldId.c_str(), fieldId.length());
    if (!field) {
        throw std::runtime_error("Xios: Null pointer for field '" + fieldId + "'");
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
    bool exists;
    cxios_field_valid_id(&exists, fieldId.c_str(), fieldId.length());
    if (exists) {
        throw std::runtime_error("Xios: Field '" + fieldId + "' already exists");
    }
    xios::CField* field = NULL;
    cxios_xml_tree_add_field(getFieldGroup(), &field, fieldId.c_str(), fieldId.length());
    if (!field) {
        throw std::runtime_error("Xios: Null pointer for field '" + fieldId + "'");
    }
    cxios_field_valid_id(&exists, fieldId.c_str(), fieldId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create field '" + fieldId + "'");
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
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_name(field)) {
        Logged::warning("Xios: Overwriting name for field '" + fieldId + "'");
    }
    cxios_set_field_name(field, fieldName.c_str(), fieldName.length());
    if (!cxios_is_defined_field_name(field)) {
        throw std::runtime_error("Xios: Failed to set name for field '" + fieldId + "'");
    }
}

/*!
 * Set the operation for a field with a given ID
 *
 * @param the field ID
 * @param operation to set
 */
void Xios::setFieldOperation(const std::string fieldId, const std::string operation)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_operation(field)) {
        Logged::warning("Xios: Overwriting operation for field '" + fieldId + "'");
    }
    cxios_set_field_operation(field, operation.c_str(), operation.length());
    if (!cxios_is_defined_field_operation(field)) {
        throw std::runtime_error("Xios: Failed to set operation for field '" + fieldId + "'");
    }
}

/*!
 * Set the grid reference for a field with a given ID
 *
 * @param the field ID
 * @param grid reference to set
 */
void Xios::setFieldGridRef(const std::string fieldId, const std::string gridRef)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_grid_ref(field)) {
        Logged::warning("Xios: Overwriting grid reference for field '" + fieldId + "'");
    }
    cxios_set_field_grid_ref(field, gridRef.c_str(), gridRef.length());
    if (!cxios_is_defined_field_grid_ref(field)) {
        throw std::runtime_error("Xios: Failed to set grid reference for field '" + fieldId + "'");
    }
}

/*!
 * Set the read access for a field with a given ID
 *
 * @param the field ID
 * @param read access to set
 */
void Xios::setFieldReadAccess(const std::string fieldId, const bool readAccess)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_read_access(field)) {
        Logged::warning("Xios: Overwriting read access for field '" + fieldId + "'");
    }
    cxios_set_field_read_access(field, readAccess);
    if (!cxios_is_defined_field_read_access(field)) {
        throw std::runtime_error("Xios: Failed to set read access for field '" + fieldId + "'");
    }
}

/*!
 * Set the frequency offset for a field with a given ID
 *
 * @param the field ID
 * @param frequency offset to set
 */
void Xios::setFieldFreqOffset(const std::string fieldId, const Duration freqOffset)
{
    xios::CField* field = getField(fieldId);
    if (cxios_is_defined_field_freq_offset(field)) {
        Logged::warning("Xios: Overwriting frequency offset for field '" + fieldId + "'");
    }
    cxios_set_field_freq_offset(field, convertDurationToXios(freqOffset));
    if (!cxios_is_defined_field_freq_offset(field)) {
        throw std::runtime_error(
            "Xios: Failed to set frequency offset for field '" + fieldId + "'");
    }
}

/*!
 * Get the name of a field with a given ID
 *
 * @param the field ID
 * @return name of the corresponding field
 */
std::string Xios::getFieldName(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_name(field)) {
        throw std::runtime_error("Xios: Undefined name for field '" + fieldId + "'");
    }
    char cStr[cStrLen];
    cxios_get_field_name(field, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the operation associated with a field with a given ID
 *
 * @param the field ID
 * @return operation used for the corresponding field
 */
std::string Xios::getFieldOperation(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_operation(field)) {
        throw std::runtime_error("Xios: Undefined operation for field '" + fieldId + "'");
    }
    char cStr[cStrLen];
    cxios_get_field_operation(field, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the grid reference associated with a field with a given ID
 *
 * @param the field ID
 * @return grid reference used for the corresponding field
 */
std::string Xios::getFieldGridRef(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_grid_ref(field)) {
        throw std::runtime_error("Xios: Undefined grid reference for field '" + fieldId + "'");
    }
    char cStr[cStrLen];
    cxios_get_field_grid_ref(field, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the read access associated with a field with a given ID
 *
 * @param the field ID
 * @return read access used for the corresponding field
 */
bool Xios::getFieldReadAccess(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_read_access(field)) {
        throw std::runtime_error("Xios: Undefined read access for field '" + fieldId + "'");
    }
    bool readAccess;
    cxios_get_field_read_access(field, &readAccess);
    return readAccess;
}

/*!
 * Get the frequency offset associated with a field with a given ID
 *
 * @param the field ID
 * @return frequency offset used for the corresponding field
 */
Duration Xios::getFieldFreqOffset(const std::string fieldId)
{
    xios::CField* field = getField(fieldId);
    if (!cxios_is_defined_field_freq_offset(field)) {
        throw std::runtime_error("Xios: Undefined frequency offset for field '" + fieldId + "'");
    }
    cxios_duration duration;
    cxios_get_field_freq_offset(field, &duration);
    char cStr[cStrLen];
    cxios_duration_convert_to_string(duration, cStr, cStrLen);
    return convertDurationFromXios(duration);
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
        throw std::runtime_error("Xios: Null pointer for group 'file_definition'");
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
    bool exists;
    cxios_file_valid_id(&exists, fileId.c_str(), fileId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Undefined file '" + fileId + "'");
    }
    xios::CFile* file = NULL;
    cxios_file_handle_create(&file, fileId.c_str(), fileId.length());
    if (!file) {
        throw std::runtime_error("Xios: Null pointer for file '" + fileId + "'");
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
    bool exists;
    cxios_file_valid_id(&exists, fileId.c_str(), fileId.length());
    if (exists) {
        throw std::runtime_error("Xios: File '" + fileId + "' already exists");
    }
    cxios_xml_tree_add_file(getFileGroup(), &file, fileId.c_str(), fileId.length());
    if (!file) {
        throw std::runtime_error("Xios: Null pointer for file '" + fileId + "'");
    }
    cxios_file_valid_id(&exists, fileId.c_str(), fileId.length());
    if (!exists) {
        throw std::runtime_error("Xios: Failed to create file '" + fileId + "'");
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
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_name(file)) {
        Logged::warning("Xios: Overwriting name for file '" + fileId + "'");
    }
    cxios_set_file_name(file, fileName.c_str(), fileName.length());
    if (!cxios_is_defined_file_name(file)) {
        throw std::runtime_error("Xios: Failed to set name for file '" + fileId + "'");
    }
}

/*!
 * Set the type of a file with a given ID
 *
 * @param the file ID
 * @param file type to set
 */
void Xios::setFileType(const std::string fileId, const std::string fileType)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_type(file)) {
        Logged::warning("Xios: Overwriting type for file '" + fileId + "'");
    }
    cxios_set_file_type(file, fileType.c_str(), fileType.length());
    if (!cxios_is_defined_file_type(file)) {
        throw std::runtime_error("Xios: Failed to set type for file '" + fileId + "'");
    }
}

/*!
 * Set the output frequency of a file with a given ID
 *
 * @param the file ID
 * @param output frequency to set
 */
void Xios::setFileOutputFreq(const std::string fileId, const Duration freq)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_output_freq(file)) {
        Logged::warning("Xios: Overwriting output frequency for file '" + fileId + "'");
    }
    cxios_set_file_output_freq(file, convertDurationToXios(freq));
    if (!cxios_is_defined_file_output_freq(file)) {
        throw std::runtime_error("Xios: Failed to set output frequency for file '" + fileId + "'");
    }
}

/*!
 * Set the split frequency of a file with a given ID
 *
 * @param the file ID
 * @param split frequency to set
 */
void Xios::setFileSplitFreq(const std::string fileId, const Duration freq)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_split_freq(file)) {
        Logged::warning("Xios: Split frequency already set for file '" + fileId + "'");
    }
    cxios_set_file_split_freq(file, convertDurationToXios(freq));
    if (!cxios_is_defined_file_split_freq(file)) {
        throw std::runtime_error("Xios: Failed to set split frequency for file '" + fileId + "'");
    }
}

/*!
 * Set the mode of a file with a given ID
 *
 * @param the file ID
 * @param file mode to set
 */
void Xios::setFileMode(const std::string fileId, const std::string mode)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_mode(file)) {
        Logged::warning("Xios: Overwriting mode for file '" + fileId + "'");
    }
    cxios_set_file_mode(file, mode.c_str(), mode.length());
    if (!cxios_is_defined_file_mode(file)) {
        throw std::runtime_error("Xios: Failed to set mode for file '" + fileId + "'");
    }
}

/*!
 * Set the parallel access mode of a file with a given ID
 *
 * @param the file ID
 * @param parallel access mode to set
 */
void Xios::setFileParAccess(const std::string fileId, const std::string parAccess)
{
    xios::CFile* file = getFile(fileId);
    if (cxios_is_defined_file_par_access(file)) {
        Logged::warning("Xios: Overwriting parallel access for file '" + fileId + "'");
    }
    cxios_set_file_par_access(file, parAccess.c_str(), parAccess.length());
    if (!cxios_is_defined_file_par_access(file)) {
        throw std::runtime_error("Xios: Failed to set parallel access for file '" + fileId + "'");
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
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_name(file)) {
        throw std::runtime_error("Xios: Undefined name for file '" + fileId + "'");
    }
    char cStr[cStrLen];
    cxios_get_file_name(file, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the type of a file with a given ID
 *
 * @param the file ID
 * @return type of the corresponding file
 */
std::string Xios::getFileType(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_type(file)) {
        throw std::runtime_error("Xios: Undefined type for file '" + fileId + "'");
    }
    char cStr[cStrLen];
    cxios_get_file_type(file, cStr, cStrLen);
    return convertCStrToCppStr(cStr, cStrLen);
}

/*!
 * Get the output frequency of a file with a given ID
 *
 * @param the file ID
 * @return the corresponding output frequency
 */
Duration Xios::getFileOutputFreq(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_output_freq(file)) {
        throw std::runtime_error("Xios: Undefined output frequency for file '" + fileId + "'");
    }
    cxios_duration duration;
    cxios_get_file_output_freq(file, &duration);
    return convertDurationFromXios(duration);
}

/*!
 * Get the split frequency of a file with a given ID
 *
 * @param the file ID
 * @return split frequency of the corresponding file
 */
Duration Xios::getFileSplitFreq(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_split_freq(file)) {
        throw std::runtime_error("Xios: Undefined split frequency for file '" + fileId + "'");
    }
    cxios_duration duration;
    cxios_get_file_split_freq(file, &duration);
    return convertDurationFromXios(duration);
}

/*!
 * Get the mode of a file with a given ID
 *
 * @param the file ID
 * @return mode of the corresponding file
 */
std::string Xios::getFileMode(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_mode(file)) {
        throw std::runtime_error("Xios: Undefined mode for file '" + fileId + "'");
    }
    char cStr[cStrLen];
    cxios_get_file_mode(file, cStr, cStrLen);
    std::string mode(cStr, cStrLen);
    boost::algorithm::trim_right(mode);
    return mode;
}

/*!
 * Get the parallel access mode of a file with a given ID
 *
 * @param the file ID
 * @return parallel access mode of the corresponding file
 */
std::string Xios::getFileParAccess(const std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    if (!cxios_is_defined_file_par_access(file)) {
        throw std::runtime_error("Xios: Undefined parallel access for file '" + fileId + "'");
    }
    char cStr[cStrLen];
    cxios_get_file_par_access(file, cStr, cStrLen);
    std::string parAccess(cStr, cStrLen);
    boost::algorithm::trim_right(parAccess);
    return parAccess;
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

/*!
 * Receive field from XIOS server that has been read from file.
 *
 * @param field name
 * @param reference to the ModelArray containing the data to be written
 */
void Xios::read(const std::string fieldId, ModelArray& modelarray)
{
    auto ndim = modelarray.nDimensions();
    auto dims = modelarray.dimensions();
    if (ndim == 2) {
        cxios_read_data_k82(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1]);
    } else if (ndim == 3) {
        cxios_read_data_k83(
            fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0], dims[1], dims[2]);
    } else if (ndim == 4) {
        cxios_read_data_k84(fieldId.c_str(), fieldId.length(), modelarray.getData(), dims[0],
            dims[1], dims[2], dims[3]);
    } else {
        throw std::invalid_argument("Only ModelArrays of dimension 2, 3, or 4 are supported");
    }
}
}

#endif
