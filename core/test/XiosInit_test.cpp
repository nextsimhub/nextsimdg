/*!
 * @file XiosInit_test.cpp
 * @brief Initialisation Testing for XIOS
 * @date Feb 28, 2023
 * @author Dr Tom Meltzer <tdm39@cam.ac.uk>
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 *
 * Test file intended to validate the initialization of the server process of
 * the XIOS class, including syncronising the state of the Nextsim grid. This
 * file serves as a set of unit tests, the system tests are elsewhere.
 */

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/doctest.h>

#include <mpi.h>
#include "include/Xios.hpp"
#include <iostream>

#include "include/Configurator.hpp"

Nextsim::Xios *xios_handler;

int main( int argc, char* argv[] ) {

  doctest::Context context;

  // Enable xios in the 'config'
  Nextsim::Configurator::clearStreams();
  std::stringstream config;
  config << "[xios]" << std::endl
    << "enable = true" << std::endl;
  std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
  Nextsim::Configurator::addStream(std::move(pcstream));

  // (Temporary) MPI Init - Required by XIOS
  MPI_Init(&argc,&argv);

  // XIOS Class Init --- Initialises server. Unknown error if initialised per test.
  // TODO: Investigate and find workaround
  xios_handler = new Nextsim::Xios;

  int result = context.run();

  if(context.shouldExit())
    return result;

  // global clean-up
  delete xios_handler;
  MPI_Finalize();
  return result;
}

/*!
 * TestXiosInitialization
 *
 * This test should pass when the Xios server process has been initialized.
 *
 */
TEST_CASE("TestXiosInitialization")
{
  REQUIRE(xios_handler->isInitialized());
}

/*!
 * #TestXiosDefaultConfiguration
 *
 * Test that default xml config has been read successfully
 *
 */
TEST_CASE("TestXiosDefaultConfiguration")
{
  std::string datetime;

  cxios_date start = xios_handler->getCalendarStart();
  //convert cxios_date to string for comparison
  datetime  = xios_handler->convertXiosDatetimeToString(start);
  REQUIRE(datetime == "2023-03-17T17:11:00Z");

  cxios_date origin = xios_handler->getCalendarOrigin();
  //convert cxios_date to string for comparison
  datetime  = xios_handler->convertXiosDatetimeToString(origin);
  REQUIRE(datetime == "2020-01-23T00:08:15Z");

  //check all elements of cxios_duration struct 
  cxios_duration duration = xios_handler->getCalendarTimestep();
  REQUIRE(duration.year == doctest::Approx(0.0));
  REQUIRE(duration.month == doctest::Approx(0.0));
  REQUIRE(duration.day == doctest::Approx(0.0));
  REQUIRE(duration.hour == doctest::Approx(1.5));
  REQUIRE(duration.minute == doctest::Approx(0.0));
  REQUIRE(duration.second == doctest::Approx(0.0));
  REQUIRE(duration.timestep == doctest::Approx(0.0));
}

/*!
 * #TestXiosDefaultConfiguration
 *
 * Test that setters work
 *
 */
TEST_CASE("TestXiosSetters")
{
  std::string datetime;

  //get Calendar start date and modify it
  cxios_date start = xios_handler->getCalendarStart();
  start.minute = 37;
  //set new start date
  xios_handler->setCalendarStart(start);
  //get Calendar modified start date
  start = xios_handler->getCalendarStart();
  //convert cxios_date to string for comparison
  datetime  = xios_handler->convertXiosDatetimeToString(start);
  REQUIRE(datetime == "2023-03-17T17:37:00Z");

  //same steps as calendar start date
  cxios_date origin = xios_handler->getCalendarOrigin();
  origin.second = 1;
  xios_handler->setCalendarOrigin(origin);
  origin = xios_handler->getCalendarOrigin();
  datetime  = xios_handler->convertXiosDatetimeToString(origin);
  REQUIRE(datetime == "2020-01-23T00:08:01Z");

  //similar approach for timestep
  cxios_duration duration = xios_handler->getCalendarTimestep();
  duration.year = 0.5;
  xios_handler->setCalendarTimestep(duration);
  duration = xios_handler->getCalendarTimestep();
  REQUIRE(duration.year == doctest::Approx(0.5));
  REQUIRE(duration.month == doctest::Approx(0.0));
  REQUIRE(duration.day == doctest::Approx(0.0));
  REQUIRE(duration.hour == doctest::Approx(1.5));
  REQUIRE(duration.minute == doctest::Approx(0.0));
  REQUIRE(duration.second == doctest::Approx(0.0));
  REQUIRE(duration.timestep == doctest::Approx(0.0));
}
