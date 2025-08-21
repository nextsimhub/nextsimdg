/*!
 * @author Tim Spain <timothy.spain@nersc.no>
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
#include "include/ModelMPI.hpp"
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
#include <ncGroup.h>
#include <ncVar.h>

#include <filesystem>
#include <sstream>

const std::string test_files_dir = TEST_FILES_DIR;
#ifdef USE_MPI
const std::string partition_filename = test_files_dir + "/partition_metadata_2.nc";
#endif

namespace Nextsim {

TEST_SUITE_BEGIN("ConfigOutput");
#ifdef USE_MPI
MPI_TEST_CASE("Test periodic output", 2)
#else
TEST_CASE("Test periodic output")
#endif
{
    size_t nx = 10;
    size_t ny = 9;

#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance(test_comm);
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
    config << "period = 3600" << std::endl; // Output every hour
    config << "start = 2020-01-11T00:00:00Z" << std::endl; // start after 10 days
    config << "field_names = " << hiceName << "," << ciceName << "," << tsurfName << ","
           << "top_melt" << std::endl;
    config << "filename = diag%m%d.nc" << std::endl;
    config << "file_period = 86400" << std::endl; // Files every day

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

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

    ModelComponent::getStore().registerArray(Protected::H_ICE, &hice);
    ModelComponent::getStore().registerArray(Protected::C_ICE, &cice);
    ModelComponent::getStore().registerArray(Protected::H_SNOW, &hsnow);
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

    IDiagnosticOutput& ido = Module::getImplementation<IDiagnosticOutput>();
    tryConfigure(ido);

    for (size_t j = 0; j < localNY - 2 * HALOWIDTH; ++j) {
        for (size_t i = 0; i < localNX - 2 * HALOWIDTH; ++i) {
            hice(i + 1, j + 1) = 0 + 0.01 * (j * nx + (i + offsetX));
            cice(i + 1, j + 1) = 0.1 + 0.01 * (j * nx + (i + offsetX));
            hsnow(i + 1, j + 1) = 0.2 + 0.01 * (j * nx + (i + offsetX));
            tsurf(i + 1, j + 1) = 0.4 + 0.01 * (j * nx + (i + offsetX));
            topMelt(i + 1, j + 1) = 0.6 + 0.01 * (j * nx + (i + offsetX));
        }
    }
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
            meta.incrementTime(Duration(3600.));
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

    const std::string specFile = diagFiles[5];
    std::set<std::string> fields = { "hice", "cice", "tsurf", "top_melt" };

    // Read the netCDF file directly
    netCDF::NcFile ncFile(specFile, netCDF::NcFile::read);
    netCDF::NcGroup metaGroup(ncFile.getGroup(IStructure::metadataNodeName()));
    netCDF::NcGroup dataGroup(ncFile.getGroup(IStructure::dataNodeName()));

    // Read the time axis
    netCDF::NcDim timeDim = dataGroup.getDim(timeName);
    // Read the time variable
    netCDF::NcVar timeVar = dataGroup.getVar(timeName);
    REQUIRE(timeDim.getSize() == hr_day);

    std::multimap<std::string, netCDF::NcVar> vars(dataGroup.getVars());
    REQUIRE(vars.size() == fields.size() + 1 + 4); // +1 for the time variable + 4 for the coords
    for (auto field : fields) {
        REQUIRE(vars.count(field) == 1);
    }
    REQUIRE(vars.count("time") == 1);
    REQUIRE(vars.count("hsnow") == 0);

    ncFile.close();

    // Clean the testing files
    for (auto fileName : diagFiles) {
        std::filesystem::remove(fileName);
    }

    Finalizer::finalize();
}
TEST_SUITE_END();
}
