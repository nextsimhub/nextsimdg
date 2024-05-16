/*!
 * @file    Xios.cpp
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @date    Fri 23 Feb 13:43:16 GMT 2024
 * @brief   XIOS interface implementation
 * @details
 *
 * Implementation of XIOS interface
 *
 * This C++ interface is designed to implement core functionality of XIOS so
 * that it can be used in nextsimdg. It is by no means meant to be
 * feature-complete. Initially the goal is to generate most of the XIOS
 * configuration in the xml definition file `iodef.xml`. As required we will
 * add more features to the C++ interface.
 *
 * To enable XIOS in nextsim add the following lines to the config file.
 *   [xios]
 *   enable = true
 *
 */
#include <boost/date_time/posix_time/time_parsers.hpp>
#if USE_XIOS

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
 *
 */
Xios::Xios() { configure(); }

//! Destructor
Xios::~Xios() { finalize(); }

//! Finalize XIOS context once xml config has been read and calendar settings updated
void Xios::context_finalize()
{
    if (isEnabled) {
        cxios_context_close_definition();
    }
}

//! close context and finialize server
void Xios::finalize()
{
    if (isEnabled) {
        cxios_context_finalize();
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
    // check if xios is enabled in the nextsim configuration
    istringstream(Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()))
        >> std::boolalpha >> isEnabled;
    if (isEnabled) {
        configureServer();
    }
}

//! Configure calendar settings
void Xios::configureServer()
{
    // initialize XIOS Server process and store MPI communicator
    clientId = "client";
    nullComm_F = MPI_Comm_c2f(MPI_COMM_NULL);
    cxios_init_client(clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F);

    // initialize nextsim context
    contextId = "nextsim";
    cxios_context_initialize(contextId.c_str(), contextId.length(), &clientComm_F);

    // initialize nextsim calendar wrapper
    cxios_get_current_calendar_wrapper(&clientCalendar);

    // initialize rank and size
    clientComm = MPI_Comm_f2c(clientComm_F);
    MPI_Comm_rank(clientComm, &rank);
    MPI_Comm_size(clientComm, &size);
}

/*!
 * verify xios server is initialized
 *
 * @return true when xios server is initialized
 */
bool Xios::isInitialized()
{
    bool init { false };
    cxios_context_is_initialized(contextId.c_str(), contextId.length(), &init);
    return init;
}

/*!
 * return datetime as std::string using ISO 8601 format (default)
 * if `isoFormat` is true format will be  2023-03-03T17:11:00Z
 * if `isoFormat` is false format will be 2023-03-03 17:11:00
 *
 * @param datetime
 * @param isoFormat as bool
 * @return datetime as a string
 */
std::string Xios::convertXiosDatetimeToString(cxios_date datetime, bool isoFormat)
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
 * helpful utility function to print cxios date.
 *
 * @param date
 */
void Xios::printCXiosDate(cxios_date date)
{
    std::cout << " year     " << date.year << std::endl;
    std::cout << " month    " << date.month << std::endl;
    std::cout << " day      " << date.day << std::endl;
    std::cout << " hour     " << date.hour << std::endl;
    std::cout << " minute   " << date.minute << std::endl;
    std::cout << " second   " << date.second << std::endl;
}

/*!
 * helpful utility function to print cxios duration.
 *
 * @param duration
 */
void Xios::printCXiosDuration(cxios_duration duration)
{
    std::cout << " year     " << duration.year << std::endl;
    std::cout << " month    " << duration.month << std::endl;
    std::cout << " day      " << duration.day << std::endl;
    std::cout << " hour     " << duration.hour << std::endl;
    std::cout << " minute   " << duration.minute << std::endl;
    std::cout << " second   " << duration.second << std::endl;
    std::cout << " timestep " << duration.timestep << std::endl;
}

/*!
 * get calendar origin
 *
 * @return calendar origin
 */
cxios_date Xios::getCalendarOrigin()
{
    cxios_date calendar_origin;
    cxios_get_calendar_wrapper_date_time_origin(clientCalendar, &calendar_origin);
    return calendar_origin;
}

/*!
 * get calendar start date
 *
 * @return calendar start date
 */
cxios_date Xios::getCalendarStart()
{
    cxios_date calendar_start;
    cxios_get_calendar_wrapper_date_start_date(clientCalendar, &calendar_start);
    return calendar_start;
}

/*!
 * get calendar timestep
 *
 * @return calendar timestep
 */
cxios_duration Xios::getCalendarTimestep()
{
    cxios_duration calendar_timestep;
    cxios_get_calendar_wrapper_timestep(clientCalendar, &calendar_timestep);
    return calendar_timestep;
}

/*!
 * get calendar step
 *
 * @return calendar step
 */
int Xios::getCalendarStep()
{
    int step = clientCalendar->getCalendar()->getStep();
    return step;
}

/*!
 * get current calendar date
 *
 * @return current calendar date
 */
std::string Xios::getCurrentDate(bool isoFormat)
{
    cxios_date xiosDate;
    cxios_get_current_date(&xiosDate);
    std::string strDate = convertXiosDatetimeToString(xiosDate, isoFormat);
    return strDate;
}

/*!
 * set calendar origin
 *
 * @param origin
 */
void Xios::setCalendarOrigin(cxios_date origin)
{
    cxios_set_calendar_wrapper_date_time_origin(clientCalendar, origin);
}

/*!
 * set calendar start date
 *
 * @param start date
 */
void Xios::setCalendarStart(cxios_date start)
{
    cxios_set_calendar_wrapper_date_start_date(clientCalendar, start);
}

/*!
 * set calendar timestep
 *
 * @param timestep
 */
void Xios::setCalendarTimestep(cxios_duration timestep)
{
    cxios_set_calendar_wrapper_timestep(clientCalendar, timestep);
    cxios_update_calendar_timestep(clientCalendar);
}

/*!
 * update xios calendar iteration/step number
 *
 * @param current step number
 */
void Xios::updateCalendar(int stepNumber) { cxios_update_calendar(stepNumber); }

/*!
 * Get the axis associated with a given ID
 *
 * @param the axis ID
 * @return a pointer to the XIOS CAxis object
 */
xios::CAxis* Xios::getAxis(std::string axisId)
{
    xios::CAxis* axis = NULL;
    cxios_axis_handle_create(&axis, axisId.c_str(), axisId.length());
    return axis;
}

/*!
 * Get the size of a given axis (the number of global points)
 *
 * @param the axis ID
 * @return size of the corresponding axis
 */
int Xios::getAxisSize(std::string axisId)
{
    int size;
    xios::CAxis* axis = getAxis(axisId);
    cxios_get_axis_n_glo(axis, &size);
    return size;
}

/*!
 * Get the values associated with a given axis
 *
 * @param the axis ID
 * @return the corresponding values
 */
std::vector<double> Xios::getAxisValues(std::string axisId)
{
    xios::CAxis* axis = getAxis(axisId);
    int size = getAxisSize(axisId);
    double* values = new double[size];
    cxios_get_axis_value(axis, values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
}

/*!
 * Set the local longitude size for a given domain
 *
 * @param the domain ID
 * @param the local longitude size
 */
void Xios::setDomainLongitudeSize(std::string domainId, int size)
{
    xios::CDomain* domain = getDomain(domainId);
    cxios_set_domain_ni(domain, size);
}

/*!
 * Set the local latitude size for a given domain
 *
 * @param the domain ID
 * @param the local longitude size
 */
void Xios::setDomainLatitudeSize(std::string domainId, int size)
{
    xios::CDomain* domain = getDomain(domainId);
    cxios_set_domain_nj(domain, size);
}

/*!
 * Set the local start longitude for a given domain
 *
 * @param the domain ID
 * @return the local start longitude
 */
void Xios::setDomainLongitudeStart(std::string domainId, int start)
{
    xios::CDomain* domain = getDomain(domainId);
    cxios_set_domain_ibegin(domain, start);
}

/*!
 * Set the local start latitude for a given domain
 *
 * @param the domain ID
 * @return the local start latitude
 */
void Xios::setDomainLatitudeStart(std::string domainId, int start)
{
    xios::CDomain* domain = getDomain(domainId);
    cxios_set_domain_jbegin(domain, start);
}

/*!
 * Set the local longitude values for a given domain
 *
 * @param the domain ID
 * @return the local longitude values
 */
void Xios::setDomainLongitudeValues(std::string domainId, std::vector<double> values)
{
    int size = getDomainLongitudeSize(domainId);
    xios::CDomain* domain = getDomain(domainId);
    cxios_set_domain_lonvalue_1d(domain, values.data(), &size);
}

/*!
 * Set the local latitude values for a given domain
 *
 * @param the domain ID
 * @return the local latitude values
 */
void Xios::setDomainLatitudeValues(std::string domainId, std::vector<double> values)
{
    int size = getDomainLatitudeSize(domainId);
    xios::CDomain* domain = getDomain(domainId);
    cxios_set_domain_latvalue_1d(domain, values.data(), &size);
}

/*!
 * Get the domain associated with a given ID
 *
 * @param the domain ID
 * @return a pointer to the XIOS CDomain object
 */
xios::CDomain* Xios::getDomain(std::string domainId)
{
    xios::CDomain* domain = NULL;
    cxios_domain_handle_create(&domain, domainId.c_str(), domainId.length());
    return domain;
}

/*!
 * Get the type of a given domain
 *
 * @param the domain ID
 * @return the corresponding domain type
 */
std::string Xios::getDomainType(std::string domainId)
{
    int size = 20;
    char cStr[size];
    xios::CDomain* domain = getDomain(domainId);
    cxios_get_domain_type(domain, cStr, size);
    std::string domainType(cStr, size);
    boost::algorithm::trim_right(domainType);
    return domainType;
}

/*!
 * Get the global longitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding global longitude size
 */
int Xios::getDomainGlobalLongitudeSize(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int size;
    cxios_get_domain_ni_glo(domain, &size);
    return size;
}

/*!
 * Get the global latitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding global latitude size
 */
int Xios::getDomainGlobalLatitudeSize(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int size;
    cxios_get_domain_nj_glo(domain, &size);
    return size;
}

/*!
 * Get the local longitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding local longitude size
 */
int Xios::getDomainLongitudeSize(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int size;
    cxios_get_domain_ni(domain, &size);
    return size;
}

/*!
 * Get the local latitude size for a given domain
 *
 * @param the domain ID
 * @return the corresponding local latitude size
 */
int Xios::getDomainLatitudeSize(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int size;
    cxios_get_domain_nj(domain, &size);
    return size;
}

/*!
 * Get the local starting longitude for a given domain
 *
 * @param the domain ID
 * @return the local starting longitude of the corresponding domain
 */
int Xios::getDomainLongitudeStart(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int start;
    cxios_get_domain_ibegin(domain, &start);
    return start;
}

/*!
 * Get the local starting latitude for a given domain
 *
 * @param the domain ID
 * @return the local starting latitude of the corresponding domain
 */
int Xios::getDomainLatitudeStart(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int start;
    cxios_get_domain_jbegin(domain, &start);
    return start;
}

/*!
 * Get the local longitude values for a given domain
 *
 * @param the domain ID
 * @return the local longitude values of the corresponding domain
 */
std::vector<double> Xios::getDomainLongitudeValues(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int size = getDomainLongitudeSize(domainId);
    double* values = new double[size];
    cxios_get_domain_lonvalue_1d(domain, values, &size);
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
std::vector<double> Xios::getDomainLatitudeValues(std::string domainId)
{
    xios::CDomain* domain = getDomain(domainId);
    int size = getDomainLatitudeSize(domainId);
    double* values = new double[size];
    cxios_get_domain_latvalue_1d(domain, values, &size);
    std::vector<double> vec(values, values + size);
    delete[] values;
    return vec;
}

/*!
 * Get the file associated with a given ID
 *
 * @param the file ID
 * @return a pointer to the XIOS CFile object
 */
xios::CFile* Xios::getFile(std::string fileId)
{
    xios::CFile* file = NULL;
    cxios_file_handle_create(&file, fileId.c_str(), fileId.length());
    return file;
}

/*!
 * Verify whether a given file ID is valid
 *
 * @param the file ID
 * @return `true` if the file ID is valid, otherwise `false`
 */
bool Xios::validFileId(std::string fileId)
{
    bool valid;
    cxios_file_valid_id(&valid, fileId.c_str(), fileId.length());
    return valid;
}
/*!
 * Verify whether an output frequency has been defined for a given file ID
 *
 * @param the file ID
 * @return `true` if the output frequency has been set, otherwise `false`
 */
bool Xios::isDefinedOutputFreq(std::string fileId)
{
    xios::CFile* file = getFile(fileId);
    return cxios_is_defined_file_output_freq(file);
}

/*!
 * Get the name of a file with a given ID
 *
 * @param the file ID
 * @return name of the corresponding file
 */
std::string Xios::getFileName(std::string fileId)
{
    int size = 20;
    char cStr[size];
    xios::CFile* file = getFile(fileId);
    cxios_get_file_name(file, cStr, size);
    std::string fileName(cStr, size);
    boost::algorithm::trim_right(fileName);
    return fileName;
}

/*!
 * Get the type of a file with a given ID
 *
 * @param the file ID
 * @return type of the corresponding file
 */
std::string Xios::getFileType(std::string fileId)
{
    int size = 20;
    char cStr[size];
    xios::CFile* file = getFile(fileId);
    cxios_get_file_type(file, cStr, size);
    std::string fileType(cStr, size);
    boost::algorithm::trim_right(fileType);
    return fileType;
}

/*!
 * Get the output frequency of a file with a given ID
 *
 * @param the file ID
 * @return the corresponding output frequency
 */
std::string Xios::getFileOutputFreq(std::string fileId)
{
    cxios_duration duration;
    xios::CFile* file = getFile(fileId);
    cxios_get_file_output_freq(file, &duration);
    int size = 20;
    char cStr[size];
    cxios_duration_convert_to_string(duration, cStr, size);
    std::string outputFreq(cStr, size);
    boost::algorithm::trim_right(outputFreq);
    return outputFreq;
}

/*!
 * send 2D field to xios server to be written to file.
 *
 * @param field name
 * @param data to be written
 * @param size of 1st dimension
 * @param size of 2nd dimension
 */
void Xios::write(const std::string fieldId, double* data, const int ni, const int nj)
{
    cxios_write_data_k82(fieldId.c_str(), fieldId.length(), data, ni, nj, -1);
}

/*!
 * send 3D field to xios server to be written to file.
 *
 * @param field name
 * @param data to be written
 * @param size of 1st dimension
 * @param size of 2nd dimension
 * @param size of 3rd dimension
 */
void Xios::write(const std::string fieldId, double* data, const int ni, const int nj, const int nk)
{
    cxios_write_data_k83(fieldId.c_str(), fieldId.length(), data, ni, nj, nk, -1);
}
}

#endif
