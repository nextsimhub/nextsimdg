/*!
 * @file    Xios.hpp
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @date    Fri 23 Feb 13:43:16 GMT 2024
 * @brief   XIOS interface header
 * @details
 *
 * Header file for XIOS interface
 *
 */
#ifndef SRC_INCLUDE_XIOS_HPP
#define SRC_INCLUDE_XIOS_HPP

#include "date.hpp"
#if USE_XIOS

#include "Configured.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <include/xios_c_interface.hpp>
#include <mpi.h>

namespace Nextsim {

class Xios : public Configured<Xios> {
public:
    Xios();
    ~Xios();

    void finalize();
    void context_finalize();
    bool isInitialized();

    void configure() override;
    void configureServer();
    void configureCalendar();

    /* Date and duration */
    xios::CDate getCurrentDate();
    std::string convertXiosDatetimeToString(cxios_date datetime, bool isoFormat = true);
    void printCXiosDate(cxios_date date);
    void printCXiosDuration(cxios_duration duration);

    /* Calendar */
    cxios_date getCalendarOrigin();
    cxios_date getCalendarStart();
    cxios_duration getCalendarTimestep();
    void getCalendarConfiguration();
    int getCalendarStep();
    std::string getCalendarDate(bool isoFormat = true);
    void setCalendarOrigin(cxios_date origin);
    void setCalendarStart(cxios_date start);
    void setCalendarTimestep(cxios_duration timestep);
    void updateCalendar(int stepNumber);

    /* I/O */
    void write(const std::string fieldstr, double* data, const int ni, const int nj);

    enum {
        ENABLED_KEY,
    };

    int rank { 0 };
    int size { 0 };
    xios::CCalendarWrapper* clientCalendar;

protected:
    bool isConfigured;

private:
    bool isEnabled;

    MPI_Comm clientComm;
    MPI_Fint clientComm_F;
    MPI_Fint nullComm_F;
    std::string clientId;
    std::string contextId;
};

}

#endif // USE_XIOS
#endif
