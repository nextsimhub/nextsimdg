/*!
 * @file Xios.hpp
 * @date 7th February 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 */

#ifndef SRC_INCLUDE_XIOS_HPP
#define SRC_INCLUDE_XIOS_HPP

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
      bool isInitialized();

      void configure() override;
      void configureServer();
      void configureCalendar();

      cxios_date getCalendarOrigin();
      cxios_date getCalendarStart();
      cxios_duration getCalendarTimestep();
      void setCalendarOrigin(cxios_date origin);
      void setCalendarStart(cxios_date start);
      void setCalendarTimestep(cxios_duration timestep);

      void getCalendarConfiguration();

      std::string getCalendarDate(bool isoFormat = true);

      void updateCalendar(int stepNumber);

      std::string convertXiosDatetimeToString(cxios_date datetime, bool isoFormat = true);

      void printCXiosDate(cxios_date date);
      void printCXiosDuration(cxios_duration duration);

      enum {
        ENABLED_KEY,
      };

      int rank{0};
      int size{0};

    protected:
      bool isConfigured;

    private:
      bool isEnabled;

      xios::CCalendarWrapper* clientCalendar;
      MPI_Comm clientComm;
      MPI_Fint clientComm_F;
      MPI_Fint nullComm_F;
      std::string clientId;
      std::string contextId;

  };

} /* end namespace Nextsim */

#endif // USE_XIOS
#endif
