/*!
 * @file ParaGrid_test.cpp
 *
 * @date 04 Aug 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#include "ModelArray.hpp"
#include <cstdlib>
#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#undef INFO
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/Finalizer.hpp"
#ifdef USE_MPI
#include "ModelMPI.hpp"
#include "include/Halo.hpp"
#endif
#include "include/IStructure.hpp"
#include "include/NextsimModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/ParametricGrid.hpp"
#include "include/gridNames.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <ncAtt.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

const std::string testFilesDir = TEST_FILES_DIR;
const std::string filename = testFilesDir + "/paraGrid_test.nc";
const std::string diagFile = "paraGrid_diag.nc";
const std::string dateString = "2000-01-01T00:00:00Z";
#ifdef USE_MPI
const std::string partitionFilename = testFilesDir + "/partition_metadata_2.nc";
#endif

static const int DG = 3;
static const int DGSTRESS = 6;
static const int CG = 2;

const size_t nx = 10;
const size_t ny = 9;
const double yFactor = 0.01;
const double xFactor = 0.1;
const double scale = 1e5;

namespace Nextsim {

void initializeTestData(
    HField& frac, DGField& fracDG, HField& mask, VertexField& coordinates, HField& x, HField& y)
{
    auto dimX = ModelArray::Dimension::X;
    auto startX = ModelArray::definedDimensions.at(dimX).start;
    auto localNX = ModelArray::definedDimensions.at(dimX).localLength;
    auto dimY = ModelArray::Dimension::Y;
    auto localNY = ModelArray::definedDimensions.at(dimY).localLength;

    // In the following loops we only set the "inner" data

    // Hfield's, DGField and mask
    for (size_t j = 0; j < localNY - 2 * Halo::haloWidth; ++j) {
        double yy = scale * (j - float(ny) / 2);
        for (size_t i = 0; i < localNX - 2 * Halo::haloWidth; ++i) {
            double xx = scale * ((i + startX) - float(nx) / 2);
            x(i + Halo::haloWidth, j + Halo::haloWidth) = xx;
            y(i + Halo::haloWidth, j + Halo::haloWidth) = yy;
            frac(i + Halo::haloWidth, j + Halo::haloWidth) = j * yFactor + (i + startX) * xFactor;
            mask(i + Halo::haloWidth, j + Halo::haloWidth) = i + startX > j ? 0 : 1;
            for (size_t d = 0; d < DG; ++d) {
                fracDG.components({ i + Halo::haloWidth, j + Halo::haloWidth })[d]
                    = frac(i + Halo::haloWidth, j + Halo::haloWidth) + d;
            }
        }
    }

    dimX = ModelArray::Dimension::XVERTEX;
    startX = ModelArray::definedDimensions.at(dimX).start;
    localNX = ModelArray::definedDimensions.at(dimX).localLength;
    dimY = ModelArray::Dimension::YVERTEX;
    localNY = ModelArray::definedDimensions.at(dimY).localLength;

    // Vetex coordinates
    for (size_t i = 0; i < localNX - 2 * Halo::haloWidth; ++i) {
        for (size_t j = 0; j < localNY - 2 * Halo::haloWidth; ++j) {
            double x = (i + startX) - 0.5 - float(nx) / 2;
            double y = j - 0.5 - float(ny) / 2;
            coordinates.components({ i + Halo::haloWidth, j + Halo::haloWidth })[0] = x * scale;
            coordinates.components({ i + Halo::haloWidth, j + Halo::haloWidth })[1] = y * scale;
        }
    }
};

TEST_SUITE_BEGIN("ParaGrid");
#ifdef USE_MPI
MPI_TEST_CASE("Write and read a ModelState-based ParaGrid restart file", 2)
#else
TEST_CASE("Write and read a ModelState-based ParaGrid restart file")
#endif
{
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    std::filesystem::remove(filename);

    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance(partitionFilename);

    const auto localNX = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const auto offsetX = metadata.getLocalCornerX();
    const auto localNY = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const auto offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNX, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, localNX + 1, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNY, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, localNY + 1, offsetY);
#else
    auto& metadata = ModelMetadata::getInstance();

    const auto localNX = nx;
    const size_t offsetX = 0;
    const auto localNY = ny;
    const size_t offsetY = 0;
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
#endif

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);

    HField fractional(ModelArray::Type::H);
    DGField fractionalDG(ModelArray::Type::DG);
    HField mask(ModelArray::Type::H);
    VertexField coordinates(ModelArray::Type::VERTEX);
    HField x;
    HField y;

    // Initialize (resize and set to zero) all ModelArrays
    ModelArray* arrays[] = { &fractional, &fractionalDG, &mask, &coordinates, &x, &y };
    for (auto arr : arrays) {
        arr->resize();
        *arr = 0.;
    }

    // populate model arrays with dummy data
    initializeTestData(fractional, fractionalDG, mask, coordinates, x, y);

    DGField hice = fractionalDG + 10;
    DGField cice = fractionalDG + 20;
    DGField hsnow = fractionalDG + 30;
    DGField damage = fractionalDG * 0.;
    HField sss = fractional;

    REQUIRE(coordinates.components({ 3, 8 })[0] - coordinates.components({ 2, 8 })[0] == scale);
    REQUIRE(coordinates.components({ 3, 8 })[1] - coordinates.components({ 3, 7 })[1] == scale);

    HField gridAzimuth;
    double gridAzimuth0 = 45.;
    gridAzimuth = gridAzimuth0;

    ModelState state = { {
                             { maskName, mask },
                             { hiceName, hice },
                             { ciceName, cice },
                             { hsnowName, hsnow },
                             { damageName, damage },
                         },
        {} };

    // A model state to set the coordinates in the metadata object
    ModelState coordState = { {
                                  { xName, x },
                                  { yName, y },
                                  { coordsName, coordinates },
                                  { gridAzimuthName, gridAzimuth },
                              },
        {} };

    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));
    // The coordinates are passed through the metadata object as affix
    // coordinates is the correct way to add coordinates to a ModelState
    metadata.extractCoordinates(coordState);
    metadata.affixCoordinates(state);
    grid.dumpModelState(state, filename, true);

    REQUIRE(std::filesystem::exists(std::filesystem::path(filename)));

    // Reset the array dimensions to make sure that the read function gets them correct
#ifdef USE_MPI
    ModelArray::setDimension(ModelArray::Dimension::X, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1, 1, 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, 1);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1);
#endif
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    ModelState ms = gridIn.getModelState(filename);

    REQUIRE(ms.data.size() == state.data.size());

    ModelArray& hiceRef = ms.data.at(hiceName);
    REQUIRE(hiceRef.nDimensions() == 2);
    REQUIRE(hiceRef.dimensions()[0] == localNX);
    REQUIRE(hiceRef.dimensions()[1] == localNY);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(hiceRef.nComponents() == DG);

    // Here we don't bother passing the coordinate arrays through a ModelMetadata object
    ModelArray& coordRef = ms.data.at(coordsName);
    REQUIRE(coordRef.nDimensions() == 2);
    REQUIRE(coordRef.nComponents() == 2);
    REQUIRE(coordRef.dimensions()[0] == localNX + 1);
    REQUIRE(coordRef.dimensions()[1] == localNY + 1);
    REQUIRE(coordRef.components({ 3, 8 })[0] - coordRef.components({ 2, 8 })[0] == scale);
    REQUIRE(coordRef.components({ 3, 8 })[1] - coordRef.components({ 3, 7 })[1] == scale);

    REQUIRE(ms.data.count(xName) > 0);
    ModelArray& xRef = ms.data.at(xName);
    auto testa = xRef(3, 8);
    auto testb = coordRef.components({ 3, 7 })[0];
    REQUIRE(xRef(3, 8) == coordRef.components({ 3, 7 })[0] + scale / 2);

    REQUIRE(ms.data.count(yName) > 0);
    ModelArray& yRef = ms.data.at(yName);
    REQUIRE(yRef(3, 8) == coordRef.components({ 2, 8 })[1] + scale / 2);

    REQUIRE(ms.data.count(gridAzimuthName) > 0);
    REQUIRE(ms.data.at(gridAzimuthName)(1, 1) == gridAzimuth0);
    std::filesystem::remove(filename);

    Finalizer::finalize();
}

#ifdef USE_MPI
MPI_TEST_CASE("Write a diagnostic ParaGrid file", 2)
#else
TEST_CASE("Write a diagnostic ParaGrid file")
#endif
{
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    REQUIRE(Module::getImplementation<IStructure>().structureType() == "parametric_rectangular");

    std::filesystem::remove(diagFile);

    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const auto localNX = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const auto offsetX = metadata.getLocalCornerX();
    const auto localNY = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const auto offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNX, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, localNX + 1, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNY, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, localNY + 1, offsetY);
#else
    auto& metadata = ModelMetadata::getInstance();

    const auto localNX = nx;
    const size_t offsetX = 0;
    const auto localNY = ny;
    const size_t offsetY = 0;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
#endif

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);

    HField fractional(ModelArray::Type::H);
    DGField fractionalDG(ModelArray::Type::DG);
    HField mask(ModelArray::Type::H);
    VertexField coordinates(ModelArray::Type::VERTEX);
    HField x;
    HField y;

    // Initialize (resize and set to zero) all ModelArrays
    ModelArray* arrays[] = { &fractional, &fractionalDG, &mask, &coordinates, &x, &y };
    for (auto arr : arrays) {
        arr->resize();
        *arr = 0.;
    }

    // populate model arrays with dummy data
    initializeTestData(fractional, fractionalDG, mask, coordinates, x, y);

    REQUIRE(fractional.nDimensions() == 2);

    DGField hice = fractionalDG + 10;
    DGField cice = fractionalDG + 20;

    HField gridAzimuth;
    double gridAzimuth0 = 45.;
    gridAzimuth = gridAzimuth0;

    ModelState state = { {
                             { maskName, mask },
                             { hiceName, hice },
                             { ciceName, cice },
                         },
        {} };

    // A model state to set the coordinates in the metadata object
    ModelState coordState = { {
                                  { xName, x },
                                  { yName, y },
                                  { coordsName, coordinates },
                                  { gridAzimuthName, gridAzimuth },
                              },
        {} };

    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));
    // The coordinates are passed through the metadata object as affix
    // coordinates is the correct way to add coordinates to a ModelState
    metadata.extractCoordinates(coordState);
    metadata.affixCoordinates(state);

    for (int t = 1; t < 5; ++t) {
        hice += 100;
        cice += 100;
        state = { {
                      { hiceName, hice },
                      { ciceName, cice },
                      { xName, x },
                      { yName, y },
                      { coordsName, coordinates },
                      { gridAzimuthName, gridAzimuth },
                  },
            {} };
        metadata.incrementTime(Duration(3600));

        grid.dumpModelState(state, diagFile, false);
    }
    pio->close(diagFile);

    REQUIRE(std::filesystem::exists(std::filesystem::path(diagFile)));

    // What do we have in the file?
    netCDF::NcFile ncFile(diagFile, netCDF::NcFile::read);

    REQUIRE(ncFile.getGroups().size() == 3);
    netCDF::NcGroup structGrp(ncFile.getGroup(IStructure::structureNodeName()));
    netCDF::NcGroup metaGrp(ncFile.getGroup(IStructure::metadataNodeName()));
    netCDF::NcGroup dataGrp(ncFile.getGroup(IStructure::dataNodeName()));

    std::string structureType;
    structGrp.getAtt(grid.typeNodeName()).getValues(structureType);
    REQUIRE(structureType == grid.structureType());

    // test data
    REQUIRE(dataGrp.getVarCount() == 7);
    netCDF::NcVar hiceVar = dataGrp.getVar(hiceName);
    netCDF::NcVar ciceVar = dataGrp.getVar(ciceName);
    netCDF::NcVar maskVar = dataGrp.getVar(maskName);
    netCDF::NcVar timeVar = dataGrp.getVar(timeName);

    // hice
    REQUIRE(hiceVar.getDimCount() == 4);

    // coordinates
    REQUIRE(dataGrp.getVars().count(xName) > 0);
    REQUIRE(dataGrp.getVars().count(yName) > 0);
    REQUIRE(dataGrp.getVars().count(coordsName) > 0);
    REQUIRE(dataGrp.getVars().count(gridAzimuthName) > 0);

    ncFile.close();

    std::filesystem::remove(diagFile);

    Finalizer::finalize();
}

#ifndef TEST_FILE_SOURCE
#define TEST_FILE_SOURCE "."
#endif

#ifdef USE_MPI
MPI_TEST_CASE("Test array ordering", 2)
#else
TEST_CASE("Test array ordering")
#endif
{
    std::string inputFilename = std::string(TEST_FILE_SOURCE) + "/ParaGridIO_input_test.nc";

    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    REQUIRE(Module::getImplementation<IStructure>().structureType() == "parametric_rectangular");

    double xFactor = 10;

#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance(partitionFilename);

    const auto localNX = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const auto offsetX = metadata.getLocalCornerX();
    const auto localNY = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const auto offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNX, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNY, offsetY);
#else
    auto& metadata = ModelMetadata::getInstance();

    const auto localNX = nx;
    const size_t offsetX = 0;
    const auto localNY = ny;
    const size_t offsetY = 0;
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
#endif

    HField index2d(ModelArray::Type::H);
    index2d.resize();
    index2d = 0.;
    std::string fieldName = "index2d";
    std::set<std::string> fields = { fieldName };
    TimePoint time;

    ModelState state = ParaGridIO::readForcingTimeStatic(fields, time, inputFilename);
    REQUIRE(state.data.count(fieldName) > 0);
    index2d = state.data.at(fieldName);
    REQUIRE(index2d(3 + Halo::haloWidth, 5 + Halo::haloWidth) == (offsetX + 3) * xFactor + 5.);

    Finalizer::finalize();
}

#ifdef USE_MPI
MPI_TEST_CASE("Check an exception is thrown for an invalid file name", 2)
#else
TEST_CASE("Check an exception is thrown for an invalid file name")
#endif
{
    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    ModelState state;

    // MD5 hash of the current output of $ date
    std::string longRandomFilename("a44f5cc1f7934a8ae8dd03a95308745d.nc");
#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadataIn = ModelMetadata::getInstance();
#endif
    REQUIRE_THROWS(state = gridIn.getModelState(longRandomFilename));

    Finalizer::finalize();
}

#ifdef USE_MPI
MPI_TEST_CASE("Check if a file with the old dimension names can be read", 2)
#else
TEST_CASE("Check if a file with the old dimension names can be read")
#endif
{
    std::string inputFilename = std::string(TEST_FILE_SOURCE) + "/old_names.nc";

    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    REQUIRE(Module::getImplementation<IStructure>().structureType() == "parametric_rectangular");

    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    // Reset the array dimensions to make sure that the read function gets them correct
#ifdef USE_MPI
    ModelArray::setDimension(ModelArray::Dimension::X, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::XCG, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::YCG, 1, 1, 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, 1);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::XCG, 1);
    ModelArray::setDimension(ModelArray::Dimension::YCG, 1);
#endif
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);
    ModelState ms = gridIn.getModelState(inputFilename);

#ifdef USE_MPI
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();
    auto localNX = metadata.getLocalExtentX();
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == localNX + 2 * Halo::haloWidth);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[1] == ny + 2 * Halo::haloWidth);
#else
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::H)[1] == ny);
#endif

    Finalizer::finalize();
}

TEST_SUITE_END();
}
