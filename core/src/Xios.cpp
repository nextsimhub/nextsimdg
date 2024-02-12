/*!
 * @file Xios.cpp
 * @date 2023-10-12
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 * @author Dr Tom Meltzer <tdm39@cam.ac.uk>
 * @brief Handler class for interfacing Nextsim-DG with the XIOS library
 * @details
 * The Xios class is designed to abstract away any complexity of iteracting
 * with the XIOS library. It assumes XIOS 2.0. Design details are available
 * in the Nextsim-DG wiki.
 *
 * This class focuses on initialising the XIOS server and client processes and
 * syncronising calendar data with Nextsim by taking information from Nextsim
 * and creating XIOS data structures or converting strings into the XIOS format
 * and using the C-bindings to pass the data to an XIOS context.
 *
 * XIOS_IODEF_PATH
 *
 * TODO: Add XIOS Grid/Axis Class/Data Structure for storing data
 * TODO: Add XIOS IO Process
 *
 * Users can enable XIOS by adding a line to the config file e.g.
 *   [xios]
 *   enable = true
 *
 */
#include <boost/date_time/posix_time/time_parsers.hpp>
#if USE_XIOS

#include "include/Xios.hpp"

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

Xios::Xios()
{
    configure();
    if (isEnabled) {
        configureServer();
    }
}

Xios::~Xios() { }

void Xios::context_finalize()
{
    if (isEnabled) {
        cxios_context_close_definition();
    }
}

void Xios::finalize()
{
    if (isEnabled) {
        cxios_context_finalize();
        cxios_finalize();
    }
}

void Xios::configure()
{
    // check if xios is enabled in the nextsim configuration
    istringstream(Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()))
        >> std::boolalpha >> isEnabled;
}

void Xios::configureServer()
{
    // initialize XIOS Server process
    clientId = "client";
    nullComm_F = MPI_Comm_c2f(MPI_COMM_NULL);
    cxios_init_client(clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F);

    // initialize nextsim context
    clientComm = MPI_Comm_f2c(clientComm_F);
    contextId = "nextsim";
    cxios_context_initialize(contextId.c_str(), contextId.length(), &clientComm_F);
    cxios_get_current_calendar_wrapper(&clientCalendar);

    MPI_Comm_rank(clientComm, &rank);
    MPI_Comm_size(clientComm, &size);
}

cxios_date Xios::getCalendarOrigin()
{
    cxios_date calendar_origin;
    cxios_get_calendar_wrapper_date_time_origin(clientCalendar, &calendar_origin);
    return calendar_origin;
}

cxios_date Xios::getCalendarStart()
{
    cxios_date calendar_start;
    cxios_get_calendar_wrapper_date_start_date(clientCalendar, &calendar_start);
    return calendar_start;
}

xios::CDate Xios::getCurrentDate()
{
    xios::CDate calendar_date;
    calendar_date = clientCalendar->getCalendar()->getCurrentDate();
    return calendar_date;
}

cxios_duration Xios::getCalendarTimestep()
{
    cxios_duration calendar_timestep;
    cxios_get_calendar_wrapper_timestep(clientCalendar, &calendar_timestep);
    return calendar_timestep;
}

int Xios::getCalendarStep()
{
    int step = -1;
    step = clientCalendar->getCalendar()->getStep();
    return step;
}

void Xios::setCalendarOrigin(cxios_date origin)
{
    cxios_set_calendar_wrapper_date_time_origin(clientCalendar, origin);
}

void Xios::setCalendarStart(cxios_date start)
{
    cxios_set_calendar_wrapper_date_start_date(clientCalendar, start);
}

void Xios::setCalendarTimestep(cxios_duration timestep)
{
    cxios_set_calendar_wrapper_timestep(clientCalendar, timestep);
    cxios_update_calendar_timestep(clientCalendar);
}

// Advance time by making a call into XIOS library. Wrapping this method
// to hide implementation details.
void Xios::updateCalendar(int stepNumber) { cxios_update_calendar(stepNumber); }

void Xios::write(const std::string fieldstr, double* data, const int ni, const int nj)
{
    cxios_write_data_k82(fieldstr.c_str(), fieldstr.length(), data, ni, nj, -1);
}

// return datetime as std::string using ISO 8601 format (default)
// (isoFormat=true)  2023-03-03T17:11:00Z
// (isoFormat=false) 2023-03-03 17:11:00
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

void Xios::printCXiosDate(cxios_date date)
{
    std::cout << " year     " << date.year << std::endl;
    std::cout << " month    " << date.month << std::endl;
    std::cout << " day      " << date.day << std::endl;
    std::cout << " hour     " << date.hour << std::endl;
    std::cout << " minute   " << date.minute << std::endl;
    std::cout << " second   " << date.second << std::endl;
}

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

// Validation Method
bool Xios::isInitialized()
{
    bool init { false };
    cxios_context_is_initialized(contextId.c_str(), contextId.length(), &init);
    return init;
}
}

#endif
