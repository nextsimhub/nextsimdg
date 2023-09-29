/*!
 * @file Xios.cpp
 * @date 7th February 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 * @brief
 * @details Class to handle interfacing with the XIOS library
 *
 *   Note:
 *     Calendar Properties must be set here or in iodef.xml before access
 *
 */
#if USE_XIOS

#include "include/Xios.hpp"

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <include/xios_c_interface.hpp>
#include <iostream>
#include <mpi.h>

bool DEBUG = false;

namespace Nextsim {

template <>
const std::map<int, std::string> Configured<Xios>::keyMap
    = { { Xios::ENABLED_KEY, "xios.enable" } };

Xios::Xios()
{

    m_isConfigured = false;
    Xios::configure();

    std::cout << m_isEnabled << std::endl;

    if (m_isEnabled) {
        std::cout << "INIT CONFIGURE ENTRY" << std::endl;
        Xios::configureServer();
        Xios::configureCalendar(m_timestep, m_start, m_origin);
        Xios::validateConfiguration();
        m_isConfigured = true;

        std::cout << "cxios_context_close_definition" << std::endl;
        cxios_context_close_definition();

        std::cout << "INIT CONFIGURE EXIT" << std::endl;
    }
    std::cout << "XIOS INIT EXIT" << std::endl;
}

Xios::~Xios()
{
    if (m_isConfigured) {
        Nextsim::Xios::finalise();
    }
}

// Inherited functionality from Configured to determine if XIOS is enabled in config file
void Xios::configure()
{
    // Is XIOS enabled / do we want to configure it?
    istringstream(Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()))
        >> std::boolalpha >> m_isEnabled;

    // TODO: Extract other pre-set values from the config file
    // For now use temporary values
    m_start = "2023-03-03 17:11:00";
    m_origin = "1970-01-01 00:00:00";
    m_timestep = "P0-0T1:0:0";
}

// Teardown XIOS client & server processes, clear data
void Xios::finalise()
{
    if (m_isEnabled) {
        cxios_context_finalize();
        cxios_finalize();
    }
}

// Setup XIOS Server process
// TODO: Better consider distributed clients
void Xios::configureServer()
{
    clientId = "client";
    nullComm_F = MPI_Comm_c2f(MPI_COMM_NULL);
    cxios_init_client(clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F);

    m_clientComm = MPI_Comm_f2c(clientComm_F);
    contextId = "xios_context";
    cxios_context_initialize(contextId.c_str(), contextId.length(), &clientComm_F);
}

// Xios Calendar Methods

// Setup the calendar information. Start Date and Timestep are defined in the iodef xml
// and are required by XIOS library. This method allows them to be modified.
// Date & timestep strings should match ISO 8601 using the dash-delimted format
// e.g. 2023-03-03T17:11:00Z or 2023-03-03 17:11:00
void Xios::configureCalendar(std::string timestep, std::string start, std::string origin)
{

    std::cout << "CONFIGURE CALENDAR " << timestep << " " << start << " " << origin << std::endl;
    // Prepare XIOS and set member variable to xios calendar object
    cxios_get_current_calendar_wrapper(&m_clientCalendar);

    std::cout << "SET TIMESTEP " << std::endl;
    setCalendarTimestep(timestep);

    std::cout << "IF " << std::endl;
    if (!origin.empty()) {

        std::cout << "SET ORIGIN " << std::endl;
        setCalendarOrigin(origin);
    }

    std::cout << "SET START " << std::endl;
    setCalendarStart(start);

    // Report to std out
    if (DEBUG)
        getCalendarConfiguration();
}

std::string Xios::getCalendarOrigin(bool isoFormat)
{
    cxios_date dorigin;
    cxios_get_calendar_wrapper_date_time_origin(m_clientCalendar, &dorigin);
    return convertXiosDatetimeToString(dorigin, isoFormat);
}

void Xios::setCalendarOrigin(std::string dorigin_str)
{
    cxios_date dorigin = convertStringToXiosDatetime(dorigin_str);
    cxios_set_calendar_wrapper_date_time_origin(m_clientCalendar, dorigin);
}

std::string Xios::getCalendarStart(bool isoFormat)
{
    cxios_date dstart;
    cxios_get_calendar_wrapper_date_start_date(m_clientCalendar, &dstart);
    return convertXiosDatetimeToString(dstart, isoFormat);
}

void Xios::setCalendarStart(std::string dstart_str)
{
    cxios_date dstart = convertStringToXiosDatetime(dstart_str);
    cxios_set_calendar_wrapper_date_start_date(m_clientCalendar, dstart);
}

std::string Xios::getCalendarTimestep()
{
    std::cout << "GET CALENDAR TIMESTEP ENTRY" << std::endl;
    const int str_size = 20;
    char dur_str[str_size];

    std::cout << "cxios_get_calendar_wrapper_timestep" << std::endl;
    cxios_get_calendar_wrapper_timestep(m_clientCalendar, &dtime);
    std::cout << "cxios_duration_convert_to_string" << std::endl;
    cxios_duration_convert_to_string(dtime, dur_str, str_size);

    std::cout << "if (DEBUG) printCXiosDuration(dtime);" << std::endl;
    if (DEBUG)
        printCXiosDuration(dtime);

    std::cout << "GET CALENDAR TIMESTEP EXIT" << std::endl;
    return dur_str;

    // TODO: The above bugs sometimes, I think the bottom might be valid too.
    // Need to figure out where I've messed up.

    // cxios_get_current_calendar_wrapper( &m_clientCalendar );
    // cxios_duration dtime;
    // cxios_get_calendar_wrapper_timestep( m_clientCalendar, &dtime );
    // float test2 = dtime.timestep;
    // return "test";
}

void Xios::setCalendarTimestep(std::string timestep_str)
{
    dtime = convertStringToXiosDuration(timestep_str);
    cxios_set_calendar_wrapper_timestep(m_clientCalendar, dtime);
    cxios_update_calendar_timestep(m_clientCalendar);
}

// Do we want to return a data structure of strings? Origin/Start/Timestep?
void Xios::getCalendarConfiguration()
{
    // Report calendar origin/start
    int rank(0);
    if (rank == 0) {
        std::string calendar_origin = getCalendarOrigin();
        std::string calendar_start = getCalendarStart(true);
        // TODO: For some reason this only works during initialisation
        std::string calendar_timestep = getCalendarTimestep();
        std::cout << "calendar time_step = " << calendar_timestep << std::endl;

        std::cout << "calendar start_date = " << calendar_start << std::endl;
        std::cout << "calendar time_origin = " << calendar_origin << std::endl;

        // This is the current date
        std::string calendar_date = getCalendarDate();
        std::cout << "calendar date " << calendar_date << std::endl;
    }
}

// Not sure this exists but maybe it should
std::string Xios::getCalendarDate(bool isoFormat)
{
    cxios_date date;
    cxios_get_current_date(&date);
    return convertXiosDatetimeToString(date, isoFormat);
}

// Advance time by making a call into XIOS library. Wrapping this method
// to hide implementation details.
void Xios::updateCalendar(int stepNumber) { cxios_update_calendar(stepNumber); }

// Utilities

// If I need an Xios utility class then this is a candidate
// TODO: Create the reverse operation
std::string Xios::convertXiosDatetimeToString(cxios_date datetime, bool isoFormat)
{
    const int str_size = 20;
    char datetime_str[20];
    if (isoFormat) {
        if (DEBUG)
            std::cout << "IF" << std::endl;
        // TODO: Find a better way of doing this?
        boost::format fmt = boost::format("%1%-%2%-%3%T%4%:%5%:%6%Z") % datetime.year
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.month)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.day)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.hour)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.minute)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.second);
        return fmt.str();
    } else {
        if (DEBUG)
            std::cout << "ELSE" << std::endl;
        cxios_date_convert_to_string(datetime, datetime_str, str_size);
    }
    return datetime_str;
}

boost::posix_time::ptime Xios::convertStringToDatetime(std::string datetime_str)
{
    // Pre-format and remove the extra characters
    // TODO: Consider using some formatter or different package to achieve this result
    int delimiter = datetime_str.find(datetime_str, 'T');
    int terminator = datetime_str.find(datetime_str, 'Z');

    if (delimiter != std::string::npos) {
        std::cout << "Delim " << delimiter << std::endl;
        datetime_str.replace(delimiter, 1, " ");
        std::cout << "Delim Exit" << std::endl;
    }
    if (terminator != std::string::npos) {
        std::cout << "Term " << terminator << std::endl;
        datetime_str.replace(terminator, 0, "");
        std::cout << "Term Exit" << std::endl;
    }

    // Use our own converter until we can interface with cxios_date_convert_from_string

    boost::posix_time::ptime datetime = boost::posix_time::time_from_string(datetime_str);

    auto time = datetime.time_of_day();
    auto date = datetime.date();
    return datetime;
}

// Creates a struct from a string in YYYY-MM-DD HH-MM-SS format
cxios_date Xios::convertStringToXiosDatetime(std::string datetime_str)
{
    boost::posix_time::ptime datetime = convertStringToDatetime(datetime_str);
    auto time = datetime.time_of_day();
    auto date = datetime.date();

    // Create the cxios_date object - this is just a simple struct
    cxios_date dstart;
    dstart.year = date.year();
    dstart.month = date.month();
    dstart.day = date.day();
    dstart.hour = time.hours();
    dstart.minute = time.minutes();
    dstart.second = time.seconds();

    // TODO: Leverage the XIOS library functionality when this is documented
    // cxios_date datetime = cxios_date_convert_from_string(datetime_str);
    return dstart;
}

// We're going to reading strings in the form
// PNNNN-NN-NNTNN:NN:NN
// which represents duration in Years-Months-Dates:Hours:Minutes:Seconds
// where padding with zeros is not required as specified in ISO 8601.
// These will be extracted and assigned to the cxios_duration struct but
// the timestep will be left at zero. A seperate method will be created to
// assign timestep values if needed.
//
// Note some invalid entries are permitted such as PY-DDD and T1000:0:0
cxios_duration Xios::convertStringToXiosDuration(std::string timestep_str)
{

    cxios_duration dduration;

    // Validate character P exists i.e. string is declared as duration
    int isDuration = timestep_str.find('P') == 0;
    if (!isDuration) {
        std::cerr << timestep_str << " is not a duration string" << std::endl;
    }
    // Remove the P
    std::string duration_str = timestep_str.substr(1);

    int delimiter = duration_str.find('T');
    if (delimiter != std::string::npos) {
        // Process substring of hours, minutes and seconds
        std::string time_str = duration_str.substr(delimiter + 1);
        boost::posix_time::time_duration time = boost::posix_time::duration_from_string(time_str);
        dduration.hour = time.hours();
        dduration.minute = time.minutes();
        dduration.second = time.seconds();
    } else {
        // No hours, minutes or seconds declared in timestep string
        dduration.hour = 0;
        dduration.minute = 0;
        dduration.second = 0;
    }

    // Until we shift to PY-M-D format we will need to manually handle the
    // years, months and days.

    std::string day_str = duration_str.substr(0, delimiter);
    int year_delimiter = day_str.find('-');

    if (year_delimiter != std::string::npos) {
        // If you find a delimiter, assume everything before it is the year component
        dduration.year = std::stod(day_str.substr(0, year_delimiter));
    }

    // FIXME: The logic here is incorrect. Testing will highlight the deficiencies but
    // I have assumed one too many delimiters. Rethink how to distinguish between the
    // formats (remove one step?)

    int month_delimiter = day_str.find('-', year_delimiter);
    if (month_delimiter != std::string::npos) {
        // If there was a second delimiter we are probably in PY-M-DD format

        int day_delimiter = day_str.find('-', month_delimiter);
        if (day_delimiter != std::string::npos) {
            // If you find another delimeter then everything before it is months
            // and after is days i.e. PY-M-D format
            dduration.month = std::stod(
                day_str.substr(month_delimiter + 1, day_delimiter - (month_delimiter + 1)));
            dduration.day = std::stod(day_str.substr(day_delimiter + 1));
        } else {
            // If there is no additional delimiter then this is days i.e. PY-DDD
            dduration.day = std::stod(day_str.substr(year_delimiter + 1));
        }
    }

    dduration.timestep = 0;

    // TODO: This needs testing like.............
    return dduration;
}

void Xios::printCXiosDuration(cxios_duration durationStructure)
{
    std::cout << "cxios_duration inspection: " << std::endl;
    std::cout << "year " << durationStructure.year << std::endl;
    std::cout << "month " << durationStructure.month << std::endl;
    std::cout << "day " << durationStructure.day << std::endl;
    std::cout << "hour " << durationStructure.hour << std::endl;
    std::cout << "minute " << durationStructure.minute << std::endl;
    std::cout << "second " << durationStructure.second << std::endl;
    std::cout << "timestep " << durationStructure.timestep << std::endl;
}

// Xios Grid/Field Data

// Method to update the Xios data
// Model -> XIOS
void Xios::setState() { }

// Method to access the Xios data
// XIOS -> Model
void Xios::getState() { }

// Method to write state of Xios data to file
void Xios::writeStateData()
{
    // configure();
}

// Method to update the Xios data state with data from file
void Xios::readStateData()
{
    // Unsure if we want this tbh
}

// Validation Methods

bool Xios::validateConfiguration()
{
    return validateServerConfiguration() && validateCalendarConfiguration()
        && validateAxisConfiguration();
}

bool Xios::validateServerConfiguration() { return true; }

bool Xios::validateCalendarConfiguration() { return true; }

bool Xios::validateAxisConfiguration() { return true; }

}

#endif
