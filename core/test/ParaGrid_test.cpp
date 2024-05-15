/*!
 * @file ParaGrid_test.cpp
 *
 * @date Oct 27, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

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

const std::string filename = "paraGrid_test.nc";
const std::string diagFile = "paraGrid_diag.nc";

static const int DG = 3;
static const int DGSTRESS = 6;
static const int CG = 2;

namespace Nextsim {

size_t c = 0;

TEST_SUITE_BEGIN("ParaGrid");
TEST_CASE("Write and read a ModelState-based ParaGrid restart file")
{
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    std::filesystem::remove(filename);

    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Set the dimension lengths
    size_t nx = 25;
    size_t ny = 15;
    size_t nz = 3;
    NZLevels::set(nz);
    size_t nxcg = CG * nx + 1;
    size_t nycg = CG * ny + 1;

    double yFactor = 0.01;
    double xFactor = 0.0001;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
    ModelArray::setDimension(ModelArray::Dimension::XCG, nxcg);
    ModelArray::setDimension(ModelArray::Dimension::YCG, nycg);

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);

    HField fractional(ModelArray::Type::H);
    DGField fractionalDG(ModelArray::Type::DG);
    HField mask(ModelArray::Type::H);
    fractional.resize();
    fractionalDG.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            fractional(i, j) = j * yFactor + i * xFactor;
            mask(i, j)
                = (i - nx / 2) * (i - nx / 2) + (j - ny / 2) * (j - ny / 2) > (nx * ny) ? 0 : 1;
            for (size_t d = 0; d < DG; ++d) {
                fractionalDG.components({ i, j })[d] = fractional(i, j) + d;
            }
        }
    }

    DGField hice = fractionalDG + 10;
    DGField cice = fractionalDG + 20;
    HField hsnow = fractional + 30;
    ZField tice(ModelArray::Type::Z);
    tice.resize();
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        for (size_t k = 0; k < nz; ++k) {
            tice.zIndexAndLayer(i, k) = fractional[i] + 40 + k;
        }
    }

    VertexField coordinates(ModelArray::Type::VERTEX);
    coordinates.resize();
    // Planar coordinates
    double scale = 1e5;

    // Vertex coordinates
    for (size_t i = 0; i < ModelArray::definedDimensions.at(ModelArray::Dimension::XVERTEX).length;
         ++i) {
        for (size_t j = 0;
             j < ModelArray::definedDimensions.at(ModelArray::Dimension::YVERTEX).length; ++j) {
            double x = i - 0.5 - nx / 2;
            double y = j - 0.5 - ny / 2;
            coordinates.components({ i, j })[0] = x * scale;
            coordinates.components({ i, j })[1] = y * scale;
        }
    }

    REQUIRE(coordinates.components({ 12, 13 })[0] - coordinates.components({ 11, 13 })[0] == scale);
    REQUIRE(coordinates.components({ 12, 13 })[1] - coordinates.components({ 12, 12 })[1] == scale);

    HField x;
    HField y;
    x.resize();
    y.resize();
    // Element coordinates
    for (size_t j = 0; j < ModelArray::size(ModelArray::Dimension::Y); ++j) {
        double yy = scale * (j - ny / 2);
        for (size_t i = 0; i < ModelArray::size(ModelArray::Dimension::X); ++i) {
            double xx = scale * (i - nx / 2);
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

    grid.dumpModelState(state, metadata, filename, true);

    REQUIRE(std::filesystem::exists(std::filesystem::path(filename)));

    // Reset the array dimensions to make sure that the read function gets them correct
    ModelArray::setDimension(ModelArray::Dimension::X, 1);
    ModelArray::setDimension(ModelArray::Dimension::Y, 1);
    ModelArray::setDimension(ModelArray::Dimension::Z, 1);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1);
    ModelArray::setDimension(ModelArray::Dimension::XCG, 1);
    ModelArray::setDimension(ModelArray::Dimension::YCG, 1);
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    ModelState ms = gridIn.getModelState(filename);

    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == ny);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[2] == NZLevels::get());

    REQUIRE(ms.data.size() == state.data.size());

    ModelArray& ticeRef = ms.data.at(ticeName);
    REQUIRE(ModelArray::nDimensions(ModelArray::Type::Z) == 3);
    REQUIRE(ticeRef.getType() == ModelArray::Type::Z);
    REQUIRE(ticeRef.nDimensions() == 3);
    REQUIRE(ticeRef.dimensions()[0] == nx);
    REQUIRE(ticeRef.dimensions()[1] == ny);
    REQUIRE(ticeRef.dimensions()[2] == NZLevels::get());

    ModelArray& hiceRef = ms.data.at(hiceName);
    REQUIRE(hiceRef.nDimensions() == 2);
    REQUIRE(hiceRef.dimensions()[0] == nx);
    REQUIRE(hiceRef.dimensions()[1] == ny);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(hiceRef.nComponents() == DG);

    REQUIRE(ticeRef(12, 14, 1) == tice(12, 14, 1));

    // Here we don't bother passing the coordinate arrays through a ModelMetadata object
    ModelArray& coordRef = ms.data.at(coordsName);
    REQUIRE(coordRef.nDimensions() == 2);
    REQUIRE(coordRef.nComponents() == 2);
    REQUIRE(coordRef.dimensions()[0] == nx + 1);
    REQUIRE(coordRef.dimensions()[1] == ny + 1);
    REQUIRE(coordRef.components({ 12, 13 })[0] - coordRef.components({ 11, 13 })[0] == scale);
    REQUIRE(coordRef.components({ 12, 13 })[1] - coordRef.components({ 12, 12 })[1] == scale);

    REQUIRE(ms.data.count(xName) > 0);
    ModelArray& xRef = ms.data.at(xName);
    REQUIRE(xRef(12, 13) == coordRef.components({ 12, 13 })[0] + scale / 2);

    REQUIRE(ms.data.count(yName) > 0);
    ModelArray& yRef = ms.data.at(yName);
    REQUIRE(yRef(12, 13) == coordRef.components({ 12, 13 })[1] + scale / 2);

    REQUIRE(ms.data.count(gridAzimuthName) > 0);
    REQUIRE(ms.data.at(gridAzimuthName)(0, 0) == gridAzimuth0);
    std::filesystem::remove(filename);
}

TEST_CASE("Write a diagnostic ParaGrid file")
{
    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    REQUIRE(Module::getImplementation<IStructure>().structureType() == "parametric_rectangular");

    std::filesystem::remove(diagFile);

    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Set the dimension lengths
    size_t nx = 30;
    size_t ny = 20;
    size_t nz = 3;
    NZLevels::set(nz);
    size_t nxcg = CG * nx + 1;
    size_t nycg = CG * ny + 1;

    double yFactor = 0.01;
    double xFactor = 0.0001;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
    ModelArray::setDimension(ModelArray::Dimension::XCG, nxcg);
    ModelArray::setDimension(ModelArray::Dimension::YCG, nycg);

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);

    HField fractional(ModelArray::Type::H);
    DGField fractionalDG(ModelArray::Type::DG);
    HField mask(ModelArray::Type::H);
    fractional.resize();
    fractionalDG.resize();
    for (size_t j = 0; j < ny; ++j) {
        for (size_t i = 0; i < nx; ++i) {
            fractional(i, j) = j * yFactor + i * xFactor;
            mask(i, j) = fractional(i, j);
            //                = (i - nx / 2) * (i - nx / 2) + (j - ny / 2) * (j - ny / 2) > (nx *
            //                ny) ? 0 : 1;
            for (size_t d = 0; d < DG; ++d) {
                fractionalDG.components({ i, j })[d] = fractional(i, j) + d;
            }
        }
    }
    double prec = 1e-9;
    REQUIRE(fractional(12, 12) - fractional(11, 12) == doctest::Approx(xFactor).epsilon(prec));
    REQUIRE(fractional(12, 12) - fractional(12, 11) == doctest::Approx(yFactor).epsilon(prec));

    REQUIRE(fractionalDG(12, 12) - fractionalDG(11, 12) == doctest::Approx(xFactor).epsilon(prec));
    REQUIRE(fractionalDG(12, 12) - fractionalDG(12, 11) == doctest::Approx(yFactor).epsilon(prec));

    DGField hice = fractionalDG + 10;
    DGField cice = fractionalDG + 20;

    VertexField coordinates(ModelArray::Type::VERTEX);
    coordinates.resize();
    // Planar coordinates
    double scale = 1e5;

    for (size_t i = 0; i < ModelArray::definedDimensions.at(ModelArray::Dimension::XVERTEX).length;
         ++i) {
        for (size_t j = 0;
             j < ModelArray::definedDimensions.at(ModelArray::Dimension::YVERTEX).length; ++j) {
            double x = i - 0.5 - nx / 2;
            double y = j - 0.5 - ny / 2;
            coordinates.components({ i, j })[0] = x * scale;
            coordinates.components({ i, j })[1] = y * scale;
        }
    }

    REQUIRE(coordinates.components({ 12, 13 })[0] - coordinates.components({ 11, 13 })[0] == scale);
    REQUIRE(coordinates.components({ 12, 13 })[1] - coordinates.components({ 12, 12 })[1] == scale);

    HField x;
    HField y;
    x.resize();
    y.resize();
    // Element coordinates
    for (size_t j = 0; j < ModelArray::size(ModelArray::Dimension::Y); ++j) {
        double yy = scale * (j - ny / 2);
        for (size_t i = 0; i < ModelArray::size(ModelArray::Dimension::X); ++i) {
            double xx = scale * (i - nx / 2);
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

    grid.dumpModelState(state, metadata, diagFile, false);

    for (int t = 1; t < 5; ++t) {
        hice += 100;
        cice += 100;
        state = { {
                      { hiceName, hice },
                      { ciceName, cice },
                  },
            {} };
        metadata.incrementTime(Duration(3600));

        grid.dumpModelState(state, metadata, diagFile, false);
    }
    pio->close(diagFile);

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
    REQUIRE(dataGrp.getVarCount() == 8);
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

#define TO_STR(s) TO_STRI(s)
#define TO_STRI(s) #s
#ifndef TEST_FILE_SOURCE
#define TEST_FILE_SOURCE .
#endif

TEST_CASE("Test array ordering")
{
    std::string inputFilename = "ParaGridIO_input_test.nc";

    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    REQUIRE(Module::getImplementation<IStructure>().structureType() == "parametric_rectangular");

    size_t nx = 9;
    size_t ny = 11;
    NZLevels::set(1);

    double xFactor = 10;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, NZLevels::get());

    HField index2d(ModelArray::Type::H);
    index2d.resize();
    std::string fieldName = "index2d";
    std::set<std::string> fields = { fieldName };
    TimePoint time;

    ModelState state = ParaGridIO::readForcingTimeStatic(
        fields, time, TO_STR(TEST_FILE_SOURCE) + std::string("/") + inputFilename);
    REQUIRE(state.data.count(fieldName) > 0);
    index2d = state.data.at(fieldName);
    REQUIRE(index2d(3, 5) == 35);
    // And that's all that's needed
}

#undef TO_STR
#undef TO_STRI

TEST_CASE("Check an exception is thrown for an invalid file name")
{
    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    ModelState state;

    // MD5 hash of the current output of $ date
    std::string longRandomFilename("a44f5cc1f7934a8ae8dd03a95308745d.nc");
    REQUIRE_THROWS(state = gridIn.getModelState(longRandomFilename));
}

TEST_CASE("Check if a file with the old dimension names can be read")
{
    std::string inputFilename = "old_names.nc";

    Module::setImplementation<IStructure>("Nextsim::ParametricGrid");

    REQUIRE(Module::getImplementation<IStructure>().structureType() == "parametric_rectangular");

    size_t nx = 1;
    size_t ny = 1;
    NZLevels::set(1);

    ParametricGrid gridIn;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    gridIn.setIO(readIO);

    // Reset the array dimensions to make sure that the read function gets them correct
    ModelArray::setDimension(ModelArray::Dimension::X, 2);
    ModelArray::setDimension(ModelArray::Dimension::Y, 2);
    ModelArray::setDimension(ModelArray::Dimension::Z, 2);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 3);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 3);
    ModelArray::setDimension(ModelArray::Dimension::XCG, 2);
    ModelArray::setDimension(ModelArray::Dimension::YCG, 2);
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);

    ModelState ms = gridIn.getModelState(filename);

    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == ny);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[2] == NZLevels::get());

}

TEST_SUITE_END();

}
