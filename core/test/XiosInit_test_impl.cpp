#include <doctest/extensions/doctest_mpi.h>
#include "include/Xios.hpp"
#include <iostream>
#include "include/Configurator.hpp"

/*!
 * TestXiosInitialization
 *
 * This test checks the Xios core implementation
 *
 */
MPI_TEST_CASE("TestXiosInitialization",2)
{

  // Enable xios in the 'config'
  Nextsim::Configurator::clearStreams();
  std::stringstream config;
  config << "[xios]" << std::endl
    << "enable = true" << std::endl;
  std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
  Nextsim::Configurator::addStream(std::move(pcstream));

  Nextsim::Xios xios_handler;
  std::string datetime;

  REQUIRE(xios_handler.isInitialized());

  cxios_date start = xios_handler.getCalendarStart();
  datetime  = xios_handler.convertXiosDatetimeToString(start);
  REQUIRE(datetime == "2023-03-17T17:11:00Z");

  cxios_date origin = xios_handler.getCalendarOrigin();
  //convert cxios_date to string for comparison
  datetime  = xios_handler.convertXiosDatetimeToString(origin);
  REQUIRE(datetime == "2020-01-23T00:08:15Z");

  //check all elements of cxios_duration struct 
  cxios_duration duration = xios_handler.getCalendarTimestep();
  REQUIRE(duration.year == doctest::Approx(0.0));
  REQUIRE(duration.month == doctest::Approx(0.0));
  REQUIRE(duration.day == doctest::Approx(0.0));
  REQUIRE(duration.hour == doctest::Approx(1.5));
  REQUIRE(duration.minute == doctest::Approx(0.0));
  REQUIRE(duration.second == doctest::Approx(0.0));
  REQUIRE(duration.timestep == doctest::Approx(0.0));

  //get Calendar start date and modify it
  start = xios_handler.getCalendarStart();
  start.minute = 37;
  //set new start date
  xios_handler.setCalendarStart(start);
  //get Calendar modified start date
  start = xios_handler.getCalendarStart();
  //convert cxios_date to string for comparison
  datetime  = xios_handler.convertXiosDatetimeToString(start);
  REQUIRE(datetime == "2023-03-17T17:37:00Z");

  //same steps as calendar start date
  origin = xios_handler.getCalendarOrigin();
  origin.second = 1;
  xios_handler.setCalendarOrigin(origin);
  origin = xios_handler.getCalendarOrigin();
  datetime  = xios_handler.convertXiosDatetimeToString(origin);
  REQUIRE(datetime == "2020-01-23T00:08:01Z");

  //similar approach for timestep
  duration = xios_handler.getCalendarTimestep();
  duration.year = 0.5;
  xios_handler.setCalendarTimestep(duration);
  duration = xios_handler.getCalendarTimestep();
  REQUIRE(duration.year == doctest::Approx(0.5));
  REQUIRE(duration.month == doctest::Approx(0.0));
  REQUIRE(duration.day == doctest::Approx(0.0));
  REQUIRE(duration.hour == doctest::Approx(1.5));
  REQUIRE(duration.minute == doctest::Approx(0.0));
  REQUIRE(duration.second == doctest::Approx(0.0));
  REQUIRE(duration.timestep == doctest::Approx(0.0));

  xios_handler.finalize();

}
