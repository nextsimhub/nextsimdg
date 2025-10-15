/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#undef INFO
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif

#include "DiagnosticOutputModule/include/ConfigOutput.hpp"

#include "include/FileCallbackCloser.hpp"
#include "include/Finalizer.hpp"
#include "include/IStructure.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"
#include "include/NextsimModule.hpp"
#include "include/gridNames.hpp"
#ifdef USE_MPI
#include "ModelMPI.hpp"
#include "include/Halo.hpp"
#endif

#include <ncDim.h>
#include <ncFile.h>
#include <ncVar.h>

#include <filesystem>
#include <sstream>

const std::string test_files_dir = TEST_FILES_DIR;
#ifdef USE_MPI
const std::string partition_filename = test_files_dir + "/paragrid_test_partition_metadata_2.nc";
#endif

namespace Nextsim {

#ifdef USE_MPI
void runMe(const bool snapshot, MPIComm myMPIComm)
#else
void runMe(const bool snapshot)
#endif
{
    size_t nx = 10;
    size_t ny = 9;

#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance(myMPIComm);
    auto& meta = ModelMetadata::getInstance(partition_filename);

    const auto localNX = meta.getLocalExtentX() + 2 * Halo::haloWidth;
    const auto offsetX = meta.getLocalCornerX();
    const auto localNY = meta.getLocalExtentY() + 2 * Halo::haloWidth;
    const auto offsetY = meta.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNX, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, localNX + 1, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNY, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, localNY + 1, offsetY);
#else
    auto& meta = ModelMetadata::getInstance();
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);

    auto offsetX = 0;
    auto localNX = nx;
    auto localNY = ny;
#endif

    Module::Module<IDiagnosticOutput>::setImplementation("Nextsim::ConfigOutput");
    std::stringstream config;
    config << "[ConfigOutput]" << std::endl;
    config << "period = P0-0T03:00:00" << std::endl; // Output every three hours
    config << "start = 2020-01-11T00:00:00Z" << std::endl; // start after 10 days
    config << "field_names = " << hiceName << "," << ciceName << "," << tsurfName << ","
           << "top_melt" << std::endl;
    config << "filename = diag%m%d.nc" << std::endl;
    config << "file_period = 86400" << std::endl; // Files every day
    if (snapshot)
        config << "snapshots = true" << std::endl;
    else
        config << "snapshots = false" << std::endl;

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    /* We need to set the model time step in the ModelMetadata instance, but are forced to set
     * starting time and duration as well, even though it's not used. */
    const Duration timeStep = Duration(3600.);
    auto& metadata = ModelMetadata::getInstance();
    metadata.setTimes(TimePoint("2010-01-01"), Duration("P0-24T00:00:00"), timeStep);

    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    HField hice(ModelArray::Type::H);
    HField cice(ModelArray::Type::H);
    HField hsnow(ModelArray::Type::H);
    HField tsurf(ModelArray::Type::H);

    // An internal diagnostic field, not made available through the data store
    HField topMelt(ModelArray::Type::H);

    hice.resize();
    cice.resize();
    hsnow.resize();
    tsurf.resize();
    topMelt.resize();

    hice = 0.;
    cice = 0.;
    hsnow = 0.;
    tsurf = 0.;
    topMelt = 0.;

    ModelComponent::getStore().registerArray(Shared::H_ICE_DG, &hice);
    ModelComponent::getStore().registerArray(Shared::C_ICE_DG, &cice);
    ModelComponent::getStore().registerArray(Shared::H_SNOW_DG, &hsnow);
    ModelComponent::getStore().registerArray(Protected::T_SURF, &tsurf);

    // Set up the coordinates, but use arrays filled with zeros
    HField latlonData(ModelArray::Type::H);
    latlonData = 0.;
    VertexField coordsData(ModelArray::Type::VERTEX);
    coordsData = 0.;
    ModelState modelCoordinates = { {
                                        { longitudeName, latlonData },
                                        { latitudeName, latlonData },
                                        { gridAzimuthName, latlonData },
                                        { coordsName, coordsData },
                                    },
        {} };
    meta.extractCoordinates(modelCoordinates);
    meta.setTime(TimePoint("2020-01-01T00:00:00Z"));

    auto& ido = Module::getImplementation<IDiagnosticOutput>();
    tryConfigure(ido);

#ifdef USE_MPI
    // offset indices by 1 (haloWidth) so only the "inner" data is initialized
    for (size_t j = 0; j < localNY - 2 * Halo::haloWidth; ++j) {
        for (size_t i = 0; i < localNX - 2 * Halo::haloWidth; ++i) {
            hice(i + 1, j + 1) = 0 + 0.01 * (j * nx + (i + offsetX));
            cice(i + 1, j + 1) = 0.1 + 0.01 * (j * nx + (i + offsetX));
            hsnow(i + 1, j + 1) = 0.2 + 0.01 * (j * nx + (i + offsetX));
            tsurf(i + 1, j + 1) = 0.4 + 0.01 * (j * nx + (i + offsetX));
            topMelt(i + 1, j + 1) = 0.6 + 0.01 * (j * nx + (i + offsetX));
        }
    }
#else
    for (size_t j = 0; j < localNY; ++j) {
        for (size_t i = 0; i < localNX; ++i) {
            hice(i, j) = 0 + 0.01 * (j * nx + (i + offsetX));
            cice(i, j) = 0.1 + 0.01 * (j * nx + (i + offsetX));
            hsnow(i, j) = 0.2 + 0.01 * (j * nx + (i + offsetX));
            tsurf(i, j) = 0.4 + 0.01 * (j * nx + (i + offsetX));
            topMelt(i, j) = 0.6 + 0.01 * (j * nx + (i + offsetX));
        }
    }
#endif
    std::vector<std::string> diagFiles;
    const std::string pfx = "diag01";
    const std::string sfx = ".nc";
    const size_t hr_day = 24;
    for (size_t day = 1; day <= 20; ++day) {
        if (day > 10) {
            diagFiles.push_back(pfx + std::to_string(day) + sfx);
        }
        double dayIncr = 100.;
        hice += dayIncr;
        cice += dayIncr;
        hsnow += dayIncr;
        tsurf += dayIncr;
        for (size_t hour = 0; hour < hr_day; ++hour) {
            double hourIncr = 1;
            hice += hourIncr;
            cice += hourIncr;
            hsnow += hourIncr;
            ModelState state = { { { "top_melt", topMelt } }, {} };

            ido.outputState(state);
            meta.incrementTime(timeStep);
        }
    }

    // Close and finalize the files
    for (const std::string& file : diagFiles) {
        FileCallbackCloser::close(file);
    }

    // Now test that there are 10 files, correctly named, and check that one of
    // them (diag0116.nc) contains what it should.
    for (const std::string& file : diagFiles) {
        REQUIRE(std::filesystem::exists(file));
    }
    // // No output should occur before the designated start date
    REQUIRE(!std::filesystem::exists(pfx + "10" + sfx));

    const int day = 14;
    const std::string specFile = diagFiles[day - 1 - 10]; // We only write files after day 10
    std::set<std::string> fields = { "hice", "cice", "tsurf", "top_melt" };

    // Read the netCDF file directly
    netCDF::NcFile ncFile(specFile, netCDF::NcFile::read);

    // Read the time axis
    netCDF::NcDim timeDim = ncFile.getDim(timeName);
    // Read the time variable
    netCDF::NcVar timeVar = ncFile.getVar(timeName);
    REQUIRE(timeDim.getSize() == hr_day / 3);

    std::multimap<std::string, netCDF::NcVar> vars(ncFile.getVars());
    REQUIRE(vars.size() == fields.size() + 1 + 4); // +1 for the time variable + 4 for the coords
    for (const auto& field : fields) {
        REQUIRE(vars.count(field) == 1);
    }
    REQUIRE(vars.count("time") == 1);
    REQUIRE(vars.count("hsnow") == 0);

    const netCDF::NcVar& var = ncFile.getVar("cice");
    std::vector<double> conc(nx * ny * 24 / 3);
    var.getVar(&conc[0]);

    constexpr int i = 3;
    constexpr int j = 4;
    // 100 per day, 1 per hour, 0.1 per variable, 0.01 per grid point
    const double coordComponent = 0.1 + 0.01 * (j * nx + (i + offsetX));
    double expectedValue = (100 + 24) * (day - 1) + 101 + coordComponent;
    if (!snapshot) {
        // First value in the average is three hours before the output and one day increment after
        expectedValue += (100 + 24) * (day - 2) + 123 + coordComponent;
        expectedValue += (100 + 24) * (day - 1) + coordComponent;
        expectedValue /= 3; // Average over three outputs
    }
    REQUIRE(conc[j * nx + i + offsetX] == doctest::Approx(expectedValue));

    ncFile.close();

    // Clean the testing files
    for (const auto& fileName : diagFiles) {
        std::filesystem::remove(fileName);
    }

    Finalizer::finalize();
}

TEST_SUITE_BEGIN("ConfigOutput");
#ifdef USE_MPI
MPI_TEST_CASE("Test averaged output", 2) { runMe(false, test_comm); }
#else
TEST_CASE("Test averaged output") { runMe(false); }
#endif

#ifdef USE_MPI
MPI_TEST_CASE("Test snapshot output", 2) { runMe(true, test_comm); }
#else
TEST_CASE("Test snapshot output") { runMe(true); }
#endif

TEST_SUITE_END();
}
