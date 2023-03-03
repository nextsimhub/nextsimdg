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

namespace Nextsim {

//! Class to handle interfacing with the XIOS library
class Xios : public Configured<Xios> {
public:
    Xios(int argc, char* argv[]);
    ~Xios();
    
    void initialise();//int argc, char* argv[]);
    void Finalise();
    bool validateConfiguration();
    bool validateServerConfiguration();
    bool validateCalendarConfiguration();
    bool validateAxisConfiguration();

    void configure() override; 
    void configureServer(int argc, char* argv[]);
    void configureCalendar();

    //Decide if I want these two and the best output type
    std::string getCalendarOrigin(bool isoFormat = true);
    void setCalendarOrigin();
    std::string getCalendarStart(bool isoFormat = true);
    void setCalendarStart(char *dstart_str, int str_size);
    std::string getCalendarTimestep();
    void setCalendarTimestep();
    void getCalendarConfiguration();
    static std::string convertXiosDatetimeToString(cxios_date datetime, bool isoFormat);
    //cxios_date convertStringToXiosDatetime(std::string datetime);

    static void writeState();
    //Arguments TBC
    void setState();
    void getState();
    void writeStateData();
    void readStateData();

    enum {
        ENABLED_KEY,
    };

    static void convertXiosDateStringToIsoDate(std::string& dateString);

protected:
    bool m_isConfigured = false;
private:
    std::string m_enabledStr;
    MPI_Comm m_clientComm;
    xios::CCalendarWrapper* m_clientCalendar;
};

} /* end namespace Nextsim */

#endif