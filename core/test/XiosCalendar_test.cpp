/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   Tests for XIOS calendars
 * @details
 * This test is designed to test calendar functionality of the C++ interface
 * for XIOS.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Finalizer.hpp"
#include "include/Model.hpp"
#include "include/ModelMPI.hpp"
#include "include/Xios.hpp"

namespace Nextsim {

/*!
 * TestXiosCalendar
 *
 * This function tests the calendar functionality of the C++ interface for XIOS. It
 * needs to be run with 1 rank i.e.,
 *
 * `mpirun -n 1 ./testXiosCalendar_MPI1`
 *
 */
MPI_TEST_CASE("TestXiosCalendar", 1)
{
    // Enable XIOS in the 'config' and provide parameters to configure it
    enableXios();
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T18:11:00Z" << std::endl;
    config << "time_step = P0-0T01:00:00" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMetadata instance based off a partition metadata file
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance("xios_test_partition_metadata_1.nc");

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configureTime();

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());

    // --- Tests for calendar API
    // Calendar type
    REQUIRE(xiosHandler.getCalendarType() == "Gregorian");
    // Calendar origin
    REQUIRE(xiosHandler.getCalendarOrigin().format() == "1970-01-01T00:00:00Z"); // Default
    TimePoint origin("2020-01-23T00:08:15Z");
    xiosHandler.setCalendarOrigin(origin);
    REQUIRE(origin == xiosHandler.getCalendarOrigin());
    REQUIRE(origin.format() == "2020-01-23T00:08:15Z");
    // Calendar start
    REQUIRE(xiosHandler.getCalendarStart().format() == "2023-03-17T17:11:00Z"); // Set in config
    TimePoint start("2023-03-17T17:11:00Z");
    xiosHandler.setCalendarStart(start);
    REQUIRE(start == xiosHandler.getCalendarStart());
    REQUIRE(start.format() == "2023-03-17T17:11:00Z");
    // Timestep
    REQUIRE(xiosHandler.getCalendarTimestep().seconds() == 3600.0); // Read from config

    xiosHandler.close_context_definition();

    // --- Tests for getCurrentDate method
    REQUIRE(xiosHandler.getCalendarStep() == 0);
    REQUIRE(xiosHandler.getCurrentDate().format() == "2023-03-17T17:11:00Z");

    // -- Tests that the timestep is set up correctly
    xiosHandler.setCalendarStep(1);
    REQUIRE(xiosHandler.getCalendarStep() == 1);
    REQUIRE(xiosHandler.getCurrentDate().format() == "2023-03-17T18:11:00Z");

    xiosHandler.context_finalize();
    Finalizer::finalize();
}

}
