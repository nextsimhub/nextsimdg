/*!
 * @file ConfigOutput_test.cpp
 *
 * @date 11 May 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif

#include "DiagnosticOutputModule/include/ConfigOutput.hpp"

#include "include/FileCallbackCloser.hpp"
#include "include/IStructure.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"
#include "include/Module.hpp"
#include "include/NZLevels.hpp"
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

#include <sstream>
#include <filesystem>

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
    size_t nx = 2;
    size_t ny = 5;
    size_t nz = 3;
    NZLevels::set(nz);

#ifdef USE_MPI
    if (test_rank == 0) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 1, 0);
    }
    if (test_rank == 1) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 1, 1);
    }
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get(), NZLevels::get(), 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());
#endif

    std::stringstream config;
    config << "[Modules]" << std::endl;
    config << "DiagnosticOutputModule = Nextsim::ConfigOutput" << std::endl;
    config << std::endl;
    config << "[ConfigOutput]" << std::endl;
    config << "period = 3600" << std::endl; // Output every hour
    config << "start = 2020-01-11T00:00:00Z" << std::endl; // start after 10 days
    config << "field_names = " << hiceName << "," << ciceName << "," << ticeName << std::endl;
    config << "filename = diag%m%d.nc" << std::endl;
    config << "file_period = 86400" << std::endl; // Files every day

    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    HField hice(ModelArray::Type::H);
    HField cice(ModelArray::Type::H);
    HField hsnow(ModelArray::Type::H);
    ZField tice(ModelArray::Type::Z);

    hice.resize();
    cice.resize();
    hsnow.resize();
    tice.resize();

    ModelComponent::getStore().registerArray(Protected::H_ICE, &hice);
    ModelComponent::getStore().registerArray(Protected::C_ICE, &cice);
    ModelComponent::getStore().registerArray(Protected::H_SNOW, &hsnow);
    ModelComponent::getStore().registerArray(Protected::T_ICE, &tice);

    ModelMetadata meta;
    meta.setTime(TimePoint("2020-01-01T00:00:00Z"));

#ifdef USE_MPI
    meta.setMpiMetadata(test_comm);
#endif

    IDiagnosticOutput& ido = Module::getImplementation<IDiagnosticOutput>();
    tryConfigure(ido);

    auto dimX = ModelArray::Dimension::X;
    auto startX = ModelArray::definedDimensions.at(dimX).start;
    auto localNX = ModelArray::definedDimensions.at(dimX).localLength;

    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < localNX; ++i) {
                tice(i, j, k) = 0.1 * k + 0.4 + 0.01 * (j * nx + (i + startX));
            }
        }
    }
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < localNX; ++i) {
            hice(i, j) = 0 + 0.01 * (j * nx + (i + startX));
            cice(i, j) = 0.1 + 0.01 * (j * nx + (i + startX));
            hsnow(i, j) = 0.2 + 0.01 * (j * nx + (i + startX));
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
        tice += dayIncr;
        for (size_t hour = 0; hour < hr_day; ++hour) {
            double hourIncr = 1;
            hice += hourIncr;
            cice += hourIncr;
            hsnow += hourIncr;
            ModelState state;

            ido.outputState(meta);
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
    std::set<std::string> fields = { "hice", "cice", "tice" };

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
    REQUIRE(vars.size() == fields.size() + 1); // +1 for the time variable
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
}
TEST_SUITE_END();
}
