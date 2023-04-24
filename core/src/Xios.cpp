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
#include "include/Xios.hpp"

#include <iostream>
#include <mpi.h>
#include <include/xios_c_interface.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>

#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Nextsim {

template <>
const std::map<int, std::string> Configured<Xios>::keyMap = {
    { Xios::ENABLED_KEY, "xios.enable" }
};

    Xios::Xios(int argc, char* argv[]) {
        Xios::configure();
        Xios::configureServer(argc, argv);
        Xios::configureCalendar();
        Xios::validateConfiguration();
    }

    Xios::~Xios()
    {
        if (m_isConfigured) {
            Nextsim::Xios::Finalise();
        }
        
    }

    //Setup XIOS server process, calendar and parse iodel.xml
    void Xios::initialise()//int argc, char* argv[])
    {
            std::cout << "XIOS --- INITIALISE ENTRY POINT" << std::endl; 
            int rank(0);
            MPI_Comm_rank( m_clientComm, &rank );
            std::cout << "Hello XIOS from proc " << rank << std::endl;


            if (rank == 0) std::cout << "XIOS --- INITIALISE EXIT POINT" << std::endl; 
    }

    //Teardown XIOS Server Processes
    void Xios::Finalise()
    {
            std::cout << "XIOS --- FINALISE" << std::endl;
            //cxios_context_close_definition();
            cxios_context_finalize();
            cxios_finalize();

            MPI_Finalize();
            std::cout << "XIOS --- FINALISE EXIT" << std::endl;
    }

    //This does not setup XIOS this sets up the XIOS class (requirement - Define abstract method in Configured)
    void Xios::configure()
    {
        std::cout << "XIOS --- CONFIGURE ENTRY POINT" << std::endl; 
        m_enabledStr = Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string());
        std::cout << "XIOS --- Enabled: " << Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()) << std::endl;

        m_isConfigured = true;
        std::cout << "XIOS --- CONFIGURE EXIT POINT" << std::endl; 
    }

    void Xios::configureServer(int argc, char* argv[])
    {

        std::cout << "XIOS --- CONFIGURE SERVER ENTRY POINT" << std::endl; 
        //Temporary MPI setup until syncronisation with assoc branch
        MPI_Init(&argc, &argv);
        int n_ranks;
        MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
        std::cout << "MPI_RANKS " << n_ranks << std::endl;
        string clientId( "client" );

        MPI_Fint clientComm_F;
        MPI_Fint nullComm_F = MPI_Comm_c2f( MPI_COMM_NULL );

        cxios_init_client( clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F );
        m_clientComm = MPI_Comm_f2c( clientComm_F );

        string contextId( "test" );
        cxios_context_initialize( contextId.c_str(), contextId.length(), &clientComm_F );

        // int rank(0);
        // MPI_Comm_rank( m_clientComm, &rank );
        // int size(1);
        // MPI_Comm_size( m_clientComm, &size );

        std::cout << "XIOS --- CONFIGURE SERVER EXIT POINT" << std::endl; 
    }




    // Xios Calendar


    void Xios::configureCalendar() 
    {

        std::cout << "XIOS --- CONFIGURE CALENDAR ENTRY POINT" << std::endl; 
        //Set member variable to xios calendar object
        cxios_get_current_calendar_wrapper( &m_clientCalendar );

                std::cout << "Enter setCalendarTimestep" << std::endl;
        setCalendarTimestep();

        //TODO: Argument from Xios::Configure, itself from callsite in Model?
        //TODO: Take string as argument and convert to char for the extern c?
        char dstart_str[20] = "2023-03-03 17:11:00";
        std::cout << "Enter setCalendarStart" << std::endl;
        //TODO: xios::CException is thrown
        setCalendarStart(dstart_str, 20);


        //Report to std out
        getCalendarConfiguration();

        std::cout << "XIOS --- CONFIGURE CALENDAR EXIT POINT" << std::endl; 
    }

    std::string Xios::getCalendarOrigin(bool isoFormat)
    {
        cxios_date dorigin;
        //char dorigin_str[20];
        cxios_get_calendar_wrapper_date_time_origin( m_clientCalendar, &dorigin );
        //cxios_date_convert_to_string( dorigin, dorigin_str, 20);
        return convertXiosDatetimeToString(dorigin, isoFormat);
    }

    //TODO: Consider if we want a std::string or char* 
    void Xios::setCalendarOrigin(char* dstart_str, int str_size)
    {
        cxios_date dstart = convertStringToXiosDatetime(dstart_str);
        cxios_set_calendar_wrapper_date_time_origin( m_clientCalendar, dstart );
    }

    std::string Xios::getCalendarStart(bool isoFormat)
    {
        std::cout << "XIOS --- GET CALENDAR START ENTRY" << std::endl;
        cxios_date dstart;
        cxios_get_calendar_wrapper_date_start_date( m_clientCalendar, &dstart);
        return convertXiosDatetimeToString(dstart, isoFormat);
        std::cout << "XIOS --- GET CALENDAR START EXIT" << std::endl;
    }

    void Xios::setCalendarStart(char* dstart_str, int old_str_size)
    {

        std::cout << "XIOS --- SET CALENDAR START ENTRY" << std::endl;
        
        //This string must be in YYYY-MM-DD HH-MM-SS format (I think the size of each can vary but the orders and delims cant)
        //char *dstart_str_test = "2023-03-03 17-11-00";

        //int str_size = 20;//11+1; 
        //const char *dstart_str_test = new char[str_size];
        //dstart_str_test = "2023-03-03 17:11:00"; 

        //const char *dstart_str_test = "2023-03-03 17:11:00"; 
        //std::string dstart_str_test = "2023-03-03 17:11:00";
        cxios_date dstart = convertStringToXiosDatetime(dstart_str);
        cxios_set_calendar_wrapper_date_start_date( m_clientCalendar, dstart );

        //cxios_set_calendar_wrapper_date_start_date( m_clientCalendar, dstart, 20);

        std::cout << "XIOS --- SET CALENDAR START EXIT" << std::endl;
    }

    std::string Xios::getCalendarTimestep()
    {
        std::cout << "XIOS --- GET CALENDAR TIMESTEP ENTRY" << std::endl;
        const int str_size = 50;//11+1; 
        char dur_str[str_size]; 
        cxios_duration dtime;
        std::cout << "get calendar wrapper timestep " << std::endl;
        cxios_get_calendar_wrapper_timestep(m_clientCalendar, &dtime);
        std::cout << "Convert Duration To String" << std::endl; 
        cxios_duration_convert_to_string(dtime, dur_str, str_size);

        std::cout << "cxios_duration inspection: " << std::endl;
        std::cout << "year " << dtime.year << std::endl;
        std::cout << "month " << dtime.month << std::endl;
        std::cout << "day " << dtime.day << std::endl;
        std::cout << "hour " << dtime.hour << std::endl;
        std::cout << "minute " << dtime.minute << std::endl;
        std::cout << "second " << dtime.second << std::endl;
        std::cout << "timestep " << dtime.timestep << std::endl;

        std::cout << "XIOS --- GET CALENDAR TIMESTEP EXIT" << std::endl;
        return dur_str;
    }

    //Do I want arguments here or what? --- how do I link back to config argument value?
    //TODO: Consider if we want a std::string or char* 
    void Xios::setCalendarTimestep() 
    {   
        //Default timestep values
        cxios_duration dtime;
        dtime.year = 0;
        dtime.month = 0;
        dtime.day = 0;
        dtime.hour = 1;
        dtime.minute = 0;
        dtime.second = 0;
        dtime.timestep = 0;

        //__XIOS_CDuration__::fromString("1ts")

        //cxios_duration dtime = convertStringToXiosDuration(dtime_str);
        cxios_set_calendar_wrapper_timestep( m_clientCalendar, dtime );
        cxios_update_calendar_timestep( m_clientCalendar );
    }

    //Do I want an array of strings? Origin/Start/Step?
    void Xios::getCalendarConfiguration()
    {
        std::cout << "GET CALENDAR CONFIGURATION ENTRY POINT" << std::endl;
        //Report calendar origin/start
        int rank(0);
        if (rank == 0) 
        {
            std::string calendar_origin = getCalendarOrigin();
            std::string calendar_start = getCalendarStart(true);
            std::string calendar_timestep = getCalendarTimestep();
            std::cout << "calendar time_origin = " << calendar_origin << std::endl;
            std::cout << "calendar start_date = " << calendar_start << std::endl;
            std::cout << "calendar time_step = " << calendar_timestep << std::endl; 
        }

        std::cout << "GET CALENDAR CONFIGURATION EXIT" << std::endl;
    }


    // Xios Grid/Field Data


    // Method to update the Xios data
    // Model -> XIOS
    void Xios::setState()
    {

    }

    // Method to access the Xios data
    // XIOS -> Model
    void Xios::getState()
    {

    }

    // Method to write state of Xios data to file
    void Xios::writeStateData()
    {
        //configure();
    }

    // Method to update the Xios data state with data from file
    void Xios::readStateData()
    {
        //Unsure if we want this tbh
    }



    // Validation Methods

    bool Xios::validateConfiguration()
    {
        return validateServerConfiguration() && validateCalendarConfiguration() && validateAxisConfiguration();
    }

    bool Xios::validateServerConfiguration()
    {
        return true;
    }

    bool Xios::validateCalendarConfiguration()
    {
        return true;
    }

    bool Xios::validateAxisConfiguration()
    {
        return true;
    }




    // Utilities

    //If I need an Xios utility class then this is a candidate
    //TODO: Create the reverse operation
    std::string Xios::convertXiosDatetimeToString(cxios_date datetime, bool isoFormat)
    {

        std::cout << "XIOS --- CONVERT XIOS DATETIME TO STRING ENTRY" << std::endl;
        const int str_size = 20;
        char datetime_str[20];
        if (isoFormat) {
            std::cout << "IF" << std::endl;
            //TODO: Find a better way of doing this?
            //TODO: The individual elements need pre-formatting. Make a function?
            boost::format fmt = boost::format("%1%-%2%-%3%T%4%:%5%:%6%Z") 
            % datetime.year
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.month)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.day)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.hour)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.minute)
            % boost::io::group(std::setw(2), std::setfill('0'), datetime.second);
            return fmt.str();
        } else {
            std::cout << "ELSE" << std::endl;
            cxios_date_convert_to_string(datetime, datetime_str, str_size);
        }

        std::cout << "XIOS --- CONVERT XIOS DATETIME TO STRING EXIT" << std::endl;
        return datetime_str;
    }

    cxios_date Xios::convertStringToXiosDatetime(std::string datetime_str)
    {
        std::cout << "ConvertStringToXiosDatetime" << std::endl;

        //Pre-format and remove the extra characters
        //TODO: Consider using some formatter to achieve this result
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

        //Use our own converter until we can interface with cxios_date_convert_from_string

        std::cout << "boost" << std::endl;
        //boost::gregorian::date datetime(boost::gregorian::from_simple_string(datetime_str));
        boost::posix_time::ptime datetime = boost::posix_time::time_from_string(datetime_str);
        auto time = datetime.time_of_day();
        auto date = datetime.date();

        std::cout << "Boost Test" << datetime << std::endl; 
        std::cout << date.year() << " " << date.month() << " " << date.day() << " ";
        std::cout << time.hours() << " " << time.minutes() << " " << time.seconds() << std::endl;

        //Create the cxios_date object
        //TODO: Investigate if there is a suitable constructor
        cxios_date dstart;
        dstart.year = date.year();
        dstart.month = date.month();
        dstart.day = date.day();
        dstart.hour = time.hours();
        dstart.minute = time.minutes();
        dstart.second = time.seconds();

        //cxios_date datetime = cxios_date_convert_from_string(datetime_str);
        return dstart;
    }

}