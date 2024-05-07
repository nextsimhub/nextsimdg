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

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <include/xios_c_interface.hpp>
#include <iostream>
#include <mpi.h>
#include <regex>
#include <string>

namespace Nextsim {

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
void Xios::configure() { configureServer(); }

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
 * get current date
 *
 * @return current date
 */
xios::CDate Xios::getCurrentDate()
{
    xios::CDate calendar_date;
    calendar_date = clientCalendar->getCalendar()->getCurrentDate();
    return calendar_date;
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
 * send 2D field to xios server to be written to file.
 *
 * @param field name
 * @param data to be written
 * @param size of 1st dimension
 * @param size of 2nd dimension
 */
void Xios::write(const std::string fieldstr, const ModelArray& modelarray)
{
    auto dim2 = modelarray.dimensions();
    cxios_write_data_k82(
        fieldstr.c_str(), fieldstr.length(), modelarray.getData(), dim2[0], dim2[1], -1);
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
}

#endif
