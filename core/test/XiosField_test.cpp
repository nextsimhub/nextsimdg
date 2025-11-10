/*!
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   Tests for XIOS fields
 * @details
 * This test is designed to test field functionality of the C++ interface
 * for XIOS.
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Finalizer.hpp"
#include "include/Model.hpp"
#include "include/ModelMPI.hpp"
#include "include/Xios.hpp"

namespace Nextsim {

/*!
 * TestXiosField
 *
 * This function tests the field functionality of the C++ interface for XIOS. It
 * needs to be run with 3 ranks i.e.,
 *
 * `mpirun -n 3 ./testXiosField_MPI3`
 *
 */
MPI_TEST_CASE("TestXiosField", 3)
{
    std::stringstream config;
    config << "[model]" << std::endl;
    config << "start = 2023-03-17T17:11:00Z" << std::endl;
    config << "stop = 2023-03-17T18:11:00Z" << std::endl;
    config << "time_step = P0-0T01:00:00" << std::endl;
    config << "restart_file = xios_test_output.nc" << std::endl;
    config << "restart_period = P0-0T03:00:00" << std::endl;
    config << "partition_file = xios_test_partition_metadata_3.nc" << std::endl;
    config << "[XiosOutput]" << std::endl;
    config << "field_names = field_A" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Create ModelMPI instance based off the test communicator
    auto& modelMPI = ModelMPI::getInstance(test_comm);

    // Create a Model and configure it so that time options are parsed
    Model model;
    model.configureRestarts();
    model.configureTime();

    // Get the Xios singleton instance and check it's initialized
    // NOTE: The singleton is created when Xios::getInstance() is first called. In this test, this
    //       happens when the time sets set by ModelMetadata::setTime(). This occurs in the call to
    //       Model::configureTime() above.
    Xios& xiosHandler = Xios::getInstance();

    // Create an axis with two points
    xiosHandler.createAxis("axis_A");
    xiosHandler.setAxisSize("axis_A", 2);

    // Create a 1D grid comprised of the single axis
    xiosHandler.createGrid("grid_1D");
    xiosHandler.gridAddAxis("grid_1D", "axis_A");

    // --- Tests for field API
    // Field creation
    // NOTE: Fields associated with files are automatically created with the appropriate read access
    // based off the XiosInput.field_names or XiosOutput.field_names entries in the config when the
    // file is created at initialisation
    const std::string fieldId = "field_A";
    REQUIRE_THROWS_WITH(xiosHandler.createField(fieldId), "Xios: Field 'field_A' already exists");
    // Disallow creation of fields that aren't in either config section
    REQUIRE_THROWS_WITH(xiosHandler.createField("field_B"),
        "Xios: Field 'field_B' cannot be found in the XiosInput or XiosOutput config sections");
    // Grid reference
    REQUIRE_THROWS_WITH(
        xiosHandler.getFieldGridRef(fieldId), "Xios: Undefined grid reference for field 'field_A'");
    const std::string gridRef = "grid_1D";
    xiosHandler.setFieldGridRef(fieldId, gridRef);
    REQUIRE(xiosHandler.getFieldGridRef(fieldId) == gridRef);
    // Read access
    // NOTE: createFile parses the associated field names, creates corresponding fields, and calls
    // setFieldReadAccess (see above note)
    REQUIRE(!xiosHandler.getFieldReadAccess(fieldId));
    // Frequency offset
    ModelMetadata& metadata = ModelMetadata::getInstance();
    Duration freqOffset = metadata.stepLength();
    xiosHandler.setFieldFreqOffset(fieldId, freqOffset);
    // TODO: Overload == for Duration
    REQUIRE(xiosHandler.getFieldFreqOffset(fieldId).seconds() == freqOffset.seconds());

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
    Finalizer::finalize();
}
}
