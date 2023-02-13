/*!
 * @file Xios.cpp
 * @date 7th February 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 */
#include "include/Xios.hpp"

#include <iostream>
#include <mpi.h>
#include <include/xios_c_interface.hpp>

namespace Nextsim {

    void Xios::configure(int argc, char* argv[])
    {
            cout << "XIOS --- CONFIGURE ENTRY POINT" << endl; 

            MPI_Init(&argc, &argv);

            int n_ranks;
            MPI_Comm_size(MPI_COMM_WORLD, &n_ranks);
            cout << "MPI_RANKS " << n_ranks << endl;


            string clientId( "client" );
            MPI_Comm clientComm;

            MPI_Fint clientComm_F ;
            MPI_Fint nullComm_F = MPI_Comm_c2f( MPI_COMM_NULL );

            cxios_init_client( clientId.c_str(), clientId.length(), &nullComm_F, &clientComm_F );
            clientComm = MPI_Comm_f2c( clientComm_F );

            int rank(0);
            MPI_Comm_rank( clientComm, &rank );
            int size(1);
            MPI_Comm_size( clientComm, &size );

            cout << "Hello XIOS from proc " << rank << endl;

            string contextId( "test" );
            cxios_context_initialize( contextId.c_str(), contextId.length(), &clientComm_F );

            cxios_date dorigin;
            char dorigin_str[20];

            xios::CCalendarWrapper* clientCalendar;
            cxios_get_current_calendar_wrapper( &clientCalendar );
            cxios_get_calendar_wrapper_date_time_origin( clientCalendar, &dorigin );

            cxios_date_convert_to_string( dorigin, dorigin_str, 20);
            if (rank == 0) cout << "calendar time_origin = " << dorigin_str << endl;

            cxios_date dstart;
            char dstart_str[20];

            cxios_get_calendar_wrapper_date_start_date( clientCalendar, &dstart );
            
            cxios_date_convert_to_string( dstart, dstart_str, 20);
            if (rank == 0) cout << "calendar start_date = " << dstart_str << endl;

            cxios_duration dtime;
            dtime.year = 0;
            dtime.month = 0;
            dtime.day = 0;
            dtime.hour = 0;
            dtime.minute = 0;
            dtime.second = 0;
            dtime.timestep = 0;
            dtime.second = 3600;
            
            cxios_set_calendar_wrapper_timestep( clientCalendar, dtime );
            
            cxios_update_calendar_timestep( clientCalendar );

            cxios_context_finalize();
            cxios_finalize();

            MPI_Finalize();
            cout << "XIOS --- MPI FINALIZED" << endl;
    }

    void Xios::writeState()
    {
        //configure();
    }

}