/*!
 * @file Xios.cpp
 * @date 7th February 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 */
#include "include/Xios.hpp"

#include <iostream>
#include <mpi.h>
#include <include/xios_c_interface.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
//#include <boost/date_time/posix_time/posix_time.hpp>

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

    //This does not setup XIOS this sets up the XIOS class (requirement - Define abstract method in Configured)
    void Xios::configure()
    {
        cout << "XIOS --- CONFIGURE ENTRY POINT" << endl; 
        m_enabledStr = Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string());
        cout << "XIOS --- Enabled: " << Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()) << endl;

        m_isConfigured = true;
        cout << "XIOS --- CONFIGURE EXIT POINT" << endl; 
    }

    void Xios::configureServer(int argc, char* argv[])
    {

        cout << "XIOS --- CONFIGURE SERVER ENTRY POINT" << endl; 
        //Temporary MPI setup until syncronisation with assoc branch
        MPI_Init(&argc, &argv);
        int n_ranks;
        MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
        cout << "MPI_RANKS " << n_ranks << endl;
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

        cout << "XIOS --- CONFIGURE SERVER EXIT POINT" << endl; 
    }

    void Xios::configureCalendar() 
    {

        cout << "XIOS --- CONFIGURE CALENDAR ENTRY POINT" << endl; 
        //Set member variable to xios calendar object
        cxios_get_current_calendar_wrapper( &m_clientCalendar );

        //Report to std out
        //getCalendarConfiguration();

        //TODO: Argument from Xios::Configure, itself from callsite in Model?
        //TODO: Take string as argument and convert to char for the extern c?
        char dstart_str[20] = "2023-03-03 17:11:00";
        cout << "Enter setCalendarStart" << endl;
        //TODO: xios::CException is thrown
        //setCalendarStart(dstart_str, 20);
        cout << "Enter setCalendarTimestep" << endl;
        setCalendarTimestep();

        //Report to std out
        getCalendarConfiguration();

        cout << "XIOS --- CONFIGURE CALENDAR EXIT POINT" << endl; 
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
    void Xios::setCalendarOrigin()
    {
    }

    std::string Xios::getCalendarStart(bool isoFormat)
    {
        cout << "XIOS --- GET CALENDAR START ENTRY" << endl;
        cxios_date dstart;
        cxios_get_calendar_wrapper_date_start_date( m_clientCalendar, &dstart);
        return convertXiosDatetimeToString(dstart, isoFormat);
        cout << "XIOS --- GET CALENDAR START EXIT" << endl;
    }

    void Xios::setCalendarStart(char* dstart_str, int str_size)
    {

        cout << "XIOS --- SET CALENDAR START ENTRY" << endl;
        
        //This string must be in YYYY-MM-DD HH-MM-SS format (I think the size of each can vary but the orders and delims cant)
        char *dstart_str_test = "2023-03-03 17-11-00";
        cxios_date *dstart = new cxios_date;
        cout << "convert from string" << endl;
        dstart = cxios_date_convert_from_string(dstart_str_test, 20);//, dstart);
        cout << "set calendar wrapper" << endl;
        cxios_set_calendar_wrapper_date_start_date( m_clientCalendar, dstart );

        //cxios_set_calendar_wrapper_date_start_date( m_clientCalendar, dstart, 20);

        cout << "XIOS --- SET CALENDAR START EXIT" << endl;
    }

    std::string Xios::getCalendarTimestep()
    {
        cout << "XIOS --- GET CALENDAR TIMESTEP ENTRY" << endl;
        const int str_size = 50;//11+1; 
        char dur_str[str_size]; 
        cxios_duration dtime;
        cxios_get_calendar_wrapper_timestep(m_clientCalendar, dtime);
        cxios_duration_convert_to_string(dtime, dur_str, str_size);

        cout << "XIOS --- GET CALENDAR TIMESTEP EXIT" << endl;
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
        dtime.hour = 0;
        dtime.minute = 0;
        dtime.second = 0;
        dtime.timestep = 0;
        dtime.second = 3600;
        
        cxios_set_calendar_wrapper_timestep( m_clientCalendar, dtime );
        cxios_update_calendar_timestep( m_clientCalendar );
    }


    //If I need an Xios utility class then this is a candidate
    //TODO: Create the reverse operation
    std::string Xios::convertXiosDatetimeToString(cxios_date datetime, bool isoFormat)
    {

        cout << "XIOS --- CONVERT XIOS DATETIME TO STRING ENTRY" << endl;
        const int str_size = 20;
        char datetime_str[20];
        if (isoFormat) {
            cout << "IF" << endl;
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
            cout << "ELSE" << endl;
            cxios_date_convert_to_string(datetime, datetime_str, str_size);
        }

        cout << "XIOS --- CONVERT XIOS DATETIME TO STRING EXIT" << endl;
        return datetime_str;
    }

    // cxios_date Xios::convertStringToXiosDatetime(std::string datetime_str)
    // {
    //     int delimiter = std::string::find(datetime_str, 'T');
    //     int terminator = std::string::find(datetime_str. 'Z');

    //     if (delimiter) {
    //         std::string::replace
    //     }
    //     if (terminator) {
    //         std::string::replace
    //     }

    //     cxios_date datetime =cxios_date_convert_from_string
    //     return datetime;
    // }

    //Do I want an array of strings? Origin/Start/Step?
    void Xios::getCalendarConfiguration()
    {
        cout << "GET CALENDAR CONFIGURATION ENTRY POINT" << endl;
        //Report calendar origin/start
        int rank(0);
        if (rank == 0) 
        {
            std::string calendar_origin = getCalendarOrigin();
            std::string calendar_start = getCalendarStart(true);
            std::string calendar_timestep = getCalendarTimestep();
            cout << "calendar time_origin = " << calendar_origin << endl;
            cout << "calendar start_date = " << calendar_start << endl;
            cout << "calendar time_step = " << calendar_timestep << endl; 
        }

        cout << "GET CALENDAR CONFIGURATION EXIT" << endl;
    }

    //Setup XIOS server process, calendar and parse iodel.xml
    void Xios::initialise()//int argc, char* argv[])
    {
            cout << "XIOS --- INITIALISE ENTRY POINT" << endl; 
            int rank(0);
            MPI_Comm_rank( m_clientComm, &rank );
            cout << "Hello XIOS from proc " << rank << endl;


            if (rank == 0) cout << "XIOS --- INITIALISE EXIT POINT" << endl; 
    }

    void Xios::Finalise()
    {
            cout << "XIOS --- FINALISE" << endl;
            //cxios_context_close_definition();
            cxios_context_finalize();
            cxios_finalize();

            MPI_Finalize();
            cout << "XIOS --- FINALISE EXIT" << endl;
    }

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

    Xios::~Xios()
    {
        if (m_isConfigured) {
            Nextsim::Xios::Finalise();
        }
        
    }

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

}