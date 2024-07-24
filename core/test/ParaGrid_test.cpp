/*!
 * @file ParaGrid_test.cpp
 *
 * @date Oct 27, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "ModelArray.hpp"
#include <cstdlib>
#ifdef USE_MPI
#include <doctest/extensions/doctest_mpi.h>
#else
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#endif

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/NZLevels.hpp"
#include "include/ParaGridIO.hpp"
#include "include/ParametricGrid.hpp"
#include "include/StructureModule.hpp"
#include "include/gridNames.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include <ncAtt.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

const std::string test_files_dir = TEST_FILES_DIR;
const std::string filename = test_files_dir + "/paraGrid_test.nc";
const std::string diagFile = "paraGrid_diag.nc";
const std::string date_string = "2000-01-01T00:00:00Z";
#ifdef USE_MPI
const std::string partition_filename = test_files_dir + "/partition_metadata_2.nc";
#endif

static const int DG = 3;
static const int DGSTRESS = 6;
static const int CG = 2;

const size_t nx = 10;
const size_t ny = 9;
const size_t nz = 3;
const double yFactor = 0.01;
const double xFactor = 0.1;
const double scale = 1e5;

namespace Nextsim {

size_t c = 0;

void initialize_test_data(HField& hfield, DGField& dgfield, HField& mask){
    hfield.resize();
    dgfield.resize();
    mask.resize();
    auto dimX = ModelArray::Dimension::X;
    auto startX = ModelArray::definedDimensions.at(dimX).start;
    auto localNX = ModelArray::definedDimensions.at(dimX).localLength;
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < localNX; ++i) {
            hfield(i, j) = j * yFactor + (i+startX) * xFactor;
            mask(i, j) = ((i + startX) - nx / 2) * ((i + startX) - nx / 2) + (j - ny / 2) * (j - ny / 2) > (nx * ny) ? 0 : 1;
            for (size_t d = 0; d < DG; ++d) {
                dgfield.components({ i, j })[d] = hfield(i, j) + d;
            }
        }
    }
};

void initialize_test_coordinates(VertexField& coordinates){
    auto dimXVertex = ModelArray::Dimension::XVERTEX;
    auto localNXVertex = ModelArray::definedDimensions.at(dimXVertex).localLength;
    auto startXVertex = ModelArray::definedDimensions.at(dimXVertex).start;
    for (size_t i = 0; i < localNXVertex; ++i) {
        for (size_t j = 0; j < ny + 1; ++j) {
            double x = (i + startXVertex) - 0.5 - float(nx) / 2;
            double y = j - 0.5 - float(ny) / 2;
            coordinates.components({ i, j })[0] = x * scale;
            coordinates.components({ i, j })[1] = y * scale;
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

    // Set the dimension lengths
    NZLevels::set(nz);


#ifdef USE_MPI
    if (test_rank == 0) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 4, 0);
      ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, 4 + 1, 0);
    }
    if (test_rank == 1) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 6, 4);
      ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, 6 + 1, 4);
    }
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get(), NZLevels::get(), 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, ny + 1, 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
#endif

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);

    HField fractional(ModelArray::Type::H);
    DGField fractionalDG(ModelArray::Type::DG);
    HField mask(ModelArray::Type::H);
    initialize_test_data(fractional, fractionalDG, mask);

    DGField hice = fractionalDG + 10;
    DGField cice = fractionalDG + 20;
    DGField hsnow = fractionalDG + 30;
    DGField damage = fractionalDG * 0.;
    HField sss = fractional;
    ZField tice(ModelArray::Type::Z);
    tice.resize();
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        for (size_t k = 0; k < nz; ++k) {
            tice.zIndexAndLayer(i, k) = fractional[i] + 40 + k;
        }
    }

    VertexField coordinates(ModelArray::Type::VERTEX);
    initialize_test_coordinates(coordinates);

    REQUIRE(coordinates.components({ 3, 8 })[0] - coordinates.components({ 2, 8 })[0] == scale);
    REQUIRE(coordinates.components({ 3, 8 })[1] - coordinates.components({ 3, 7 })[1] == scale);

    // MPI domain is only split in x-direction for this test
    // the following will be set correctly with MPI ON and OFF
    auto dimX = ModelArray::Dimension::X;
    auto startX = ModelArray::definedDimensions.at(dimX).start;
    auto localNX = ModelArray::definedDimensions.at(dimX).localLength;
    auto dimXVertex = ModelArray::Dimension::XVERTEX;
    auto localNXVertex = ModelArray::definedDimensions.at(dimXVertex).localLength;
    auto startXVertex = ModelArray::definedDimensions.at(dimXVertex).start;

    HField x;
    HField y;
    x.resize();
    y.resize();
    // Element coordinates
    for (size_t j = 0; j < ny; ++j) {
        double yy = scale * (j - float(ny) / 2);
        for (size_t i = 0; i < localNX; ++i) {
            double xx = scale * ((i + startX) - float(nx) / 2);
            x(i, j) = xx;
            y(i, j) = yy;
        }
    }

    HField gridAzimuth;
    double gridAzimuth0 = 45.;
    gridAzimuth = gridAzimuth0;

    ModelState state = { {
                             { maskName, mask },
                             { hiceName, hice },
                             { ciceName, cice },
                             { hsnowName, hsnow },
                             { damageName, damage },
                             { ticeName, tice },
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

    ModelMetadata metadata;
    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));
    // The coordinates are passed through the metadata object as affix
    // coordinates is the correct way to add coordinates to a ModelState
    metadata.extractCoordinates(coordState);
    metadata.affixCoordinates(state);

#ifdef USE_MPI
    metadata.setMpiMetadata(test_comm);
    metadata.globalExtentX = nx;
    metadata.globalExtentY = ny;
    metadata.localCornerX = startX;
    metadata.localCornerY = 0;
    metadata.localExtentX = localNX;
    metadata.localExtentY = ny;
#endif
    grid.dumpModelState(state, metadata, filename, true);

    REQUIRE(std::filesystem::exists(std::filesystem::path(filename)));

    // Reset the array dimensions to make sure that the read function gets them correct
#ifdef USE_MPI
    ModelArray::setDimension(ModelArray::Dimension::X, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1, 1, 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, 1);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1);
    ModelArray::setDimension(ModelArray::Dimension::Z, 1);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1);
#endif
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

#ifdef USE_MPI
    ModelMetadata metadataIn(partition_filename, test_comm);
    metadataIn.setTime(TimePoint(date_string));
    ModelState ms = gridIn.getModelState(filename, metadataIn);
#else
    ModelState ms = gridIn.getModelState(filename);
#endif

    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == localNX);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == ny);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[2] == NZLevels::get());

    REQUIRE(ms.data.size() == state.data.size());

    ModelArray& ticeRef = ms.data.at(ticeName);
    REQUIRE(ModelArray::nDimensions(ModelArray::Type::Z) == 3);
    REQUIRE(ticeRef.getType() == ModelArray::Type::Z);
    REQUIRE(ticeRef.nDimensions() == 3);
    REQUIRE(ticeRef.dimensions()[0] == localNX);
    REQUIRE(ticeRef.dimensions()[1] == ny);
    REQUIRE(ticeRef.dimensions()[2] == NZLevels::get());

    ModelArray& hiceRef = ms.data.at(hiceName);
    REQUIRE(hiceRef.nDimensions() == 2);
    REQUIRE(hiceRef.dimensions()[0] == localNX);
    REQUIRE(hiceRef.dimensions()[1] == ny);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(hiceRef.nComponents() == DG);

    REQUIRE(ticeRef(4, 9, 1) == tice(4, 9, 1));

    // Here we don't bother passing the coordinate arrays through a ModelMetadata object
    ModelArray& coordRef = ms.data.at(coordsName);
    REQUIRE(coordRef.nDimensions() == 2);
    REQUIRE(coordRef.nComponents() == 2);
    REQUIRE(coordRef.dimensions()[0] == localNXVertex);
    REQUIRE(coordRef.dimensions()[1] == ny + 1);
    REQUIRE(coordRef.components({ 3, 8 })[0] - coordRef.components({ 2, 8 })[0] == scale);
    REQUIRE(coordRef.components({ 3, 8 })[1] - coordRef.components({ 3, 7 })[1] == scale);

    REQUIRE(ms.data.count(xName) > 0);
    ModelArray& xRef = ms.data.at(xName);
    REQUIRE(xRef(3, 8) == coordRef.components({ 3, 7 })[0] + scale / 2);

    REQUIRE(ms.data.count(yName) > 0);
    ModelArray& yRef = ms.data.at(yName);
    REQUIRE(yRef(3, 8) == coordRef.components({ 2, 8 })[1] + scale / 2);

    REQUIRE(ms.data.count(gridAzimuthName) > 0);
    REQUIRE(ms.data.at(gridAzimuthName)(0, 0) == gridAzimuth0);
    std::filesystem::remove(filename);
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

    NZLevels::set(nz);

#ifdef USE_MPI
    if (test_rank == 0) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 4, 0);
      ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, 4 + 1, 0);
    }
    if (test_rank == 1) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 6, 4);
      ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, 6 + 1, 4);
    }
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get(), NZLevels::get(), 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, ny + 1, 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
#endif

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);

    // MPI domain is only split in x-direction for this test
    // the following will be set correctly with MPI ON and OFF
    auto dimX = ModelArray::Dimension::X;
    auto startX = ModelArray::definedDimensions.at(dimX).start;
    auto localNX = ModelArray::definedDimensions.at(dimX).localLength;
    auto dimXVertex = ModelArray::Dimension::XVERTEX;
    auto localNXVertex = ModelArray::definedDimensions.at(dimXVertex).localLength;
    auto startXVertex = ModelArray::definedDimensions.at(dimXVertex).start;

    HField fractional(ModelArray::Type::H);
    DGField fractionalDG(ModelArray::Type::DG);
    HField mask(ModelArray::Type::H);
    initialize_test_data(fractional, fractionalDG, mask);

    REQUIRE(fractional.nDimensions() == 2);

    DGField hice = fractionalDG + 10;
    DGField cice = fractionalDG + 20;

    VertexField coordinates(ModelArray::Type::VERTEX);
    initialize_test_coordinates(coordinates);

    HField x;
    HField y;
    x.resize();
    y.resize();
    // Element coordinates
    for (size_t j = 0; j < ny; ++j) {
        double yy = scale * (j - float(ny) / 2);
        for (size_t i = 0; i < localNX; ++i) {
            double xx = scale * ((i + startX) - float(nx) / 2);
            x(i, j) = xx;
            y(i, j) = yy;
        }
    }

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


    ModelMetadata metadata;
    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));
    // The coordinates are passed through the metadata object as affix
    // coordinates is the correct way to add coordinates to a ModelState
    metadata.extractCoordinates(coordState);
    metadata.affixCoordinates(state);

#ifdef USE_MPI
    metadata.setMpiMetadata(test_comm);
    metadata.globalExtentX = nx;
    metadata.globalExtentY = ny;
    metadata.localCornerX = startX;
    metadata.localCornerY = 0;
    metadata.localExtentX = localNX;
    metadata.localExtentY = ny;
#endif

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

        grid.dumpModelState(state, metadata, diagFile, false);
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

    // TODO test metadata

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

    size_t nx = 9;
    size_t ny = 11;
    NZLevels::set(1);

    double xFactor = 10;

#ifdef USE_MPI
    if (test_rank == 0) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 4, 0);
    }
    if (test_rank == 1) {
      ModelArray::setDimension(ModelArray::Dimension::X, nx, 5, 4);
    }
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, ny, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get(), NZLevels::get(), 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());
#endif

    HField index2d(ModelArray::Type::H);
    index2d.resize();
    std::string fieldName = "index2d";
    std::set<std::string> fields = { fieldName };
    TimePoint time;

    ModelState state = ParaGridIO::readForcingTimeStatic(fields, time, inputFilename);
    REQUIRE(state.data.count(fieldName) > 0);
    index2d = state.data.at(fieldName);
    REQUIRE(index2d(3, 5) == 35);
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
    ModelMetadata metadataIn(partition_filename, test_comm);
    metadataIn.setTime(TimePoint(date_string));
    REQUIRE_THROWS(state = gridIn.getModelState(longRandomFilename, metadataIn));
#else
    REQUIRE_THROWS(state = gridIn.getModelState(longRandomFilename));
#endif

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

    size_t nx = 2;
    size_t ny = 1;
    NZLevels::set(1);

    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    // Reset the array dimensions to make sure that the read function gets them correct
#ifdef USE_MPI
    ModelArray::setDimension(ModelArray::Dimension::X, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::Z, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::XCG, 1, 1, 0);
    ModelArray::setDimension(ModelArray::Dimension::YCG, 1, 1, 0);
#else
    ModelArray::setDimension(ModelArray::Dimension::X, 1);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1);
    ModelArray::setDimension(ModelArray::Dimension::Z, 1);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::XCG, 1);
    ModelArray::setDimension(ModelArray::Dimension::YCG, 1);
#endif
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

#ifdef USE_MPI
    ModelMetadata metadata;
    metadata.setMpiMetadata(test_comm);
    if (metadata.mpiMyRank == 0) {
        metadata.localCornerX = 0;
    }
    if (metadata.mpiMyRank == 1) {
        metadata.localCornerX = 1;
    }
    metadata.globalExtentX = nx;
    metadata.globalExtentY = ny;
    metadata.localCornerY = 0;
    metadata.localExtentX = 1;
    metadata.localExtentY = ny;
    metadata.setTime(TimePoint(date_string));
    ModelState ms = gridIn.getModelState(inputFilename, metadata);
#else
    ModelState ms = gridIn.getModelState(inputFilename);
#endif

    auto localNX = ModelArray::definedDimensions.at(ModelArray::Dimension::X).localLength;
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == localNX);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == ny);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[2] == NZLevels::get());
}

TEST_SUITE_END();

}
