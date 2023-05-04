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
#include <boost/date_time/posix_time/posix_time.hpp>

namespace Nextsim {

template <>
const std::map<int, std::string> Configured<Xios>::keyMap = {
    { Xios::ENABLED_KEY, "xios.enable" }
};

    Xios::Xios(int argc, char* argv[]) {
        Xios::configure();
        if (m_isEnabled)
        {
            Xios::configureServer(argc, argv);
            Xios::configureCalendar(m_timestep, m_start, m_origin);
            Xios::validateConfiguration();
            m_isConfigured = true;
            cxios_context_close_definition();
        }
    }
    
    Xios::~Xios()
    {
        if (m_isConfigured) {
            Nextsim::Xios::finalise();
        }
        
    }

    //Inherited functionality from Configured to determine if XIOS is enabled in config file
    void Xios::configure()
    {
        //Is XIOS enabled in the config file?
        istringstream( Configured::getConfiguration(keyMap.at(ENABLED_KEY), std::string()) )
            >> std::boolalpha >> m_isEnabled;

        //TODO: Extract other pre-set values from the config file
        //For now use temporary values
        m_start = "2023-03-03 17:11:00";
        m_origin = "1970-01-01 00:00:00";
        m_timestep = "P0-0T1:0:0";
    }

    //Teardown XIOS client & server processes, clear data
    void Xios::finalise()
    {
            cxios_context_finalize();
            cxios_finalize();

            //Temporary MPI teardown until pre-requisite feature reaches develop branch
            MPI_Finalize();
    }

    // Setup XIOS Server process
    // TODO: Consider distributed clients
    void Xios::configureServer(int argc, char* argv[])
    {

        //Temporary MPI setup until pre-requisite feature reaches develop branch
        MPI_Init(&argc, &argv);

        clientId = "client";
        nullComm_F = MPI_Comm_c2f( MPI_COMM_NULL );
        cxios_init_client( clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F );

        m_clientComm = MPI_Comm_f2c( clientComm_F );
        contextId = "xios_context";
        cxios_context_initialize( contextId.c_str(), contextId.length(), &clientComm_F );
    }


    // Xios Calendar
    
    void Xios::configureCalendar(std::string timestep, std::string start, std::string origin) 
    {

        std::cout << "XIOS --- CONFIGURE CALENDAR ENTRY POINT" << std::endl; 
        //Set member variable to xios calendar object
        cxios_get_current_calendar_wrapper( &m_clientCalendar );

        std::cout << "Enter setCalendarTimestep" << std::endl;
        setCalendarTimestep(timestep);

        if (!origin.empty()) {
            setCalendarOrigin(origin);
        }

        //TODO: Argument from Xios::Configure, itself from callsite in Model?
        //TODO: Take string as argument and convert to char for the extern c?
        std::cout << "Enter setCalendarStart" << std::endl;
        //TODO: xios::CException is thrown
        setCalendarStart(start);


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
    void Xios::setCalendarOrigin(std::string dorigin_str)
    {
        std::cout << "XIOS --- SET CALENDAR ORIGIN ENTRY" << std::endl;
        //Creates a struct from a string in YYYY-MM-DD HH-MM-SS format
        cxios_date dorigin = convertStringToXiosDatetime(dorigin_str);
        cxios_set_calendar_wrapper_date_time_origin(m_clientCalendar, dorigin);
        std::cout << "XIOS --- SET CALENDAR ORIGIN EXIT" << std::endl;
    }

    std::string Xios::getCalendarStart(bool isoFormat)
    {
        std::cout << "XIOS --- GET CALENDAR START ENTRY" << std::endl;
        cxios_date dstart;
        cxios_get_calendar_wrapper_date_start_date( m_clientCalendar, &dstart);
        return convertXiosDatetimeToString(dstart, isoFormat);
        std::cout << "XIOS --- GET CALENDAR START EXIT" << std::endl;
    }

    void Xios::setCalendarStart(std::string dstart_str)
    {

        std::cout << "XIOS --- SET CALENDAR START ENTRY" << std::endl;
        
        //Creates a struct from a string in YYYY-MM-DD HH-MM-SS format
        cxios_date dstart = convertStringToXiosDatetime(dstart_str);
        cxios_set_calendar_wrapper_date_start_date(m_clientCalendar, dstart);

        std::cout << "XIOS --- SET CALENDAR START EXIT" << std::endl;
    }

    std::string Xios::getCalendarTimestep()
    {
        std::cout << "XIOS --- GET CALENDAR TIMESTEP ENTRY" << std::endl;
        
        const int str_size = 50;//11+1; 
        char dur_str[str_size]; 
        //cxios_duration dtime;
        std::cout << "get calendar wrapper timestep " << std::endl;
        cxios_get_calendar_wrapper_timestep(m_clientCalendar, &dtime);
        std::cout << "Convert Duration To String" << std::endl; 
        cxios_duration_convert_to_string(dtime, dur_str, str_size);

        printCXiosDuration(dtime);

        std::cout << "XIOS --- GET CALENDAR TIMESTEP EXIT" << std::endl;
        return dur_str;
    }

    //Do I want arguments here or what? --- how do I link back to config argument value?
    //TODO: Consider if we want a std::string or char* 
    void Xios::setCalendarTimestep(std::string timestep_str) 
    {   
        // //Default timestep values
        // cxios_duration dtime;
        // dtime.year = 0;
        // dtime.month = 0;
        // dtime.day = 0;
        // dtime.hour = 1;
        // dtime.minute = 0;
        // dtime.second = 0;
        // dtime.timestep = 0;

        dtime = convertStringToXiosDuration(timestep_str);
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
            std::cout << "calendar time_origin = " << calendar_origin << std::endl;

            std::string calendar_start = getCalendarStart(true);
            std::cout << "calendar start_date = " << calendar_start << std::endl;

            // For some reason this only works during initialisation
            //std::string calendar_timestep = getCalendarTimestep();
            //std::cout << "calendar time_step = " << calendar_timestep << std::endl; 

            std::string calendar_date = getCalendarDate();
            std::cout << "calendar date " << calendar_date << std::endl; 
        }
        std::cout << "GET CALENDAR CONFIGURATION EXIT" << std::endl;
        
    }

    //Not sure this exists but maybe it should
    std::string Xios::getCalendarDate(bool isoFormat)
    {
        std::cout << "XIOS --- GET CALENDAR DATE ENTRY" << std::endl;
        cxios_date date;
        cxios_get_current_date(&date);
        return convertXiosDatetimeToString(date, isoFormat);
    }

    // Advance time by making a call into XIOS library. Wrapping this method
    // to hide implementation details.
    void Xios::updateCalendar(int stepNumber)
    {
        cxios_update_calendar(stepNumber);
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

    boost::posix_time::ptime Xios::convertStringToDatetime(std::string datetime_str)
    {
        //Pre-format and remove the extra characters
        //The intention is to use this or something like this for duration too
        std::cout << "ConvertStringToDatetime" << std::endl;
        //TODO: Consider using some formatter or different package to achieve this result
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

        return datetime;
    }

    cxios_date Xios::convertStringToXiosDatetime(std::string datetime_str)
    {
        std::cout << "ConvertStringToXiosDatetime" << std::endl;
        boost::posix_time::ptime datetime = convertStringToDatetime(datetime_str);
        auto time = datetime.time_of_day();
        auto date = datetime.date();
        //Create the cxios_date object
        //TODO: Investigate if there is a suitable constructor
        //cxios_date dstart(date.year(), date.month(), date.day(), time.hours() ...
        //                    time.minutes(), time.seconds() );
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

    // We're going to reading strings in the form
    // PNNNN-NN-NNTNN:NN:NN 
    // which represents duration in Years-Months-Dates:Hours:Minutes:Seconds
    // where padding with zeros is not required as specified in ISO 8601.
    // These will be extracted and assigned to the cxios_duration struct but
    // the timestep will be left at zero. A seperate method will be created to
    // assign timestep values if needed. 
    //
    // Note some invalid entries are permitted such as PY-DDD and T1000:0:0
    cxios_duration Xios::convertStringToXiosDuration(std::string timestep_str) {
        std::cout << "ConvertStringToXiosDuration " << timestep_str << " " << timestep_str.find('P') << std::endl;

        //Validate P exists
        int isDuration = timestep_str.find('P') == 0;
        if (!isDuration) 
        {
            //TODO: Error
            std::cout << timestep_str << " is not a duration string" << std::endl;
        }

        cxios_duration dduration;

        
        //Remove the P, denoting duration
        std::string duration_str = timestep_str.substr(1);

        std::cout << "Duration str " << duration_str << std::endl;

        // Until we shift to PY-M-D format we will need to manually handle the
        // years, months and days.

        int delimiter = duration_str.find('T');


        std::cout << "Duration str '" << duration_str << "'. Delim: " << delimiter << std::endl;

        if (delimiter != std::string::npos) {
            // we have minutes hours and seconds
            std::string time_str = duration_str.substr(delimiter+1);
            //boost::posix_time::ptime datetime = convertStringToDatetime(time_str);
            //auto time = datetime.time_of_day();

            boost::posix_time::time_duration time = boost::posix_time::duration_from_string(time_str);

            dduration.hour = time.hours();
            dduration.minute = time.minutes();
            dduration.second = time.seconds();
        } else {
            dduration.hour = 0;
            dduration.minute = 0;
            dduration.second = 0;
        }

        std::cout << "Duration Hour: " << dduration.hour << std::endl; 
        std::cout << "Duration Minute: " << dduration.minute << std::endl; 
        std::cout << "Duration Second: " << dduration.second << std::endl; 

        std::string day_str = duration_str.substr(0,delimiter);
        int year_delimiter = day_str.find('-');

        if (year_delimiter != std::string::npos) {
            //std::cout << day_str.substr(0,year_delimiter) << std::endl;
            dduration.year = std::stod( day_str.substr(0,year_delimiter) );
        }

        int month_delimiter = day_str.find('-', year_delimiter);
        if (month_delimiter != std::string::npos) {

            int day_delimiter = day_str.find('-', month_delimiter);
            if (day_delimiter != std::string::npos) {
                //If you find another delimeter then everything before it is months
                //and after is days i.e. PY-M-D format
                dduration.month = std::stod( 
                    day_str.substr(month_delimiter+1,day_delimiter-(month_delimiter+1)) );
                dduration.day = std::stod( day_str.substr(day_delimiter+1) );
            } else {
                //If there is no additional delimiter then this is days i.e. PY-DDD
                dduration.day = std::stod( day_str.substr(year_delimiter+1) );
            }

        }

        std::cout << "Duration Year: " << dduration.year << std::endl; 
        std::cout << "Duration Month: " << dduration.month << std::endl; 
        std::cout << "Duration Day: " << dduration.day << std::endl; 

        //TODO: This needs testing like............. 
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

}