/*!
 * @file Xios.hpp
 * @date 7th February 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 */

#ifndef SRC_INCLUDE_XIOS_HPP
#define SRC_INCLUDE_XIOS_HPP

#include <mpi.h>
#include <include/xios_c_interface.hpp>
#include "Configured.hpp"
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Nextsim {

//! Class to handle interfacing with the XIOS library
class Xios : public Configured<Xios> {
public:
    Xios(int argc, char* argv[]);
    ~Xios();
    
    //void initialise();//int argc, char* argv[]);
    void finalise();
    bool validateConfiguration();
    bool validateServerConfiguration();
    bool validateCalendarConfiguration();
    bool validateAxisConfiguration();

    void configure() override; 
    void configureServer(int argc, char* argv[]);
    void configureCalendar();

    //Decide if I want these two and the best output type
    std::string getCalendarOrigin(bool isoFormat = true);
    void setCalendarOrigin(char *dstart_str, int str_size);
    std::string getCalendarStart(bool isoFormat = true);
    void setCalendarStart(char *dstart_str, int str_size);
    std::string getCalendarTimestep();
    void setCalendarTimestep(std::string timestep_str);
    
    void getCalendarConfiguration();

    std::string getCalendarDate(bool isoFormat = true);
    void updateCalendar(int stepNumber);
    
    static std::string convertXiosDatetimeToString(cxios_date datetime, bool isoFormat);
    static boost::posix_time::ptime convertStringToDatetime(std::string datetime);
    static cxios_date convertStringToXiosDatetime(std::string datetime);
    static cxios_duration convertStringToXiosDuration(std::string duration);

    void printCXiosDuration(cxios_duration durationStructure);

    static void writeState();
    //Arguments TBC
    void setState();
    void getState();
    void writeStateData();
    void readStateData();

    enum {
        ENABLED_KEY,
    };

    //TODO: Doesn't Exist? -> Remove
    static void convertXiosDateStringToIsoDate(std::string& dateString);

protected:
    bool m_isConfigured = false;
private:
    bool m_isEnabled;

    xios::CCalendarWrapper* m_clientCalendar;
    MPI_Comm m_clientComm;
    MPI_Fint clientComm_F;
    MPI_Fint nullComm_F;
    std::string clientId;
    std::string contextId;

    cxios_duration dtime;
};

} /* end namespace Nextsim */

#endif