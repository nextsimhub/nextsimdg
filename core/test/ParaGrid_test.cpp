/*!
 * @file ParaGrid_test.cpp
 *
 * @date Oct 27, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/Configurator.hpp"
#include "include/ConfiguredModule.hpp"
#include "include/ParaGridIO.hpp"
#include "include/ParametricGrid.hpp"
#include "include/gridNames.hpp"
#include "include/IStructureModule.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>

#include <ncAtt.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>

#include <iostream> // FIXME remove me

const std::string filename = "paraGrid_test.nc";
const std::string diagFile = "paraGrid_diag.nc";

static const int DG = 3;
static const int DGSTRESS = 6;
static const int CG = 2;

namespace Nextsim {

size_t c = 0;

TEST_CASE("Write and read a ModelState-based ParaGrid restart file", "[ParametricGrid]")
{
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    Module::setImplementation<IStructure>("ParametricGrid");
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    std::FILE* lun = std::fopen(filename.c_str(), "r");
    if (lun != NULL) {
        std::remove(filename.c_str());
    }
    std::fclose(lun);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ParametricGrid grid;
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ParaGridIO* pio = new ParaGridIO(grid);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    grid.setIO(pio);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    // Set the dimension lengths
    size_t nx = 25;
    size_t ny = 15;
    size_t nz = 3;
    size_t nxcg = CG * nx + 1;
    size_t nycg = CG * ny + 1;

    double yFactor = 0.01;
    double xFactor = 0.0001;
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::XCG, nxcg);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::YCG, nycg);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelArray::setNComponents(ModelArray::Type::DG, DG);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setNComponents(ModelArray::Type::DGSTRESS, DGSTRESS);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setNComponents(ModelArray::Type::VERTEX, ModelArray::nCoords);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    HField fractional(ModelArray::Type::H);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    DGField fractionalDG(ModelArray::Type::DG);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    HField mask(ModelArray::Type::H);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    fractional.resize();
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    fractionalDG.resize();
    std::cerr << "PGt checkpoint " << ++c << std::endl;
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
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    DGField hice = fractionalDG + 10;
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    DGField cice = fractionalDG + 20;
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    DGField hsnow = fractionalDG + 30;
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ZField tice(ModelArray::Type::Z);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    tice.resize();
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    for (size_t i = 0; i < ModelArray::size(ModelArray::Type::H); ++i) {
        for (size_t k = 0; k < nz; ++k) {
            tice.zIndexAndLayer(i, k) = fractional[i] + 40 + k;
        }
    }
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    VertexField coordinates(ModelArray::Type::VERTEX);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    coordinates.resize();
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    // Planar coordinates
    double scale = 1e5;
    std::cerr << "PGt checkpoint " << ++c << std::endl;

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
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    REQUIRE(coordinates.components({ 12, 13 })[0] - coordinates.components({ 11, 13 })[0] == scale);
    REQUIRE(coordinates.components({ 12, 13 })[1] - coordinates.components({ 12, 12 })[1] == scale);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelState state = { {
                             { maskName, mask },
                             { hiceName, hice },
                             { ciceName, cice },
                             { hsnowName, hsnow },
                             { ticeName, tice },
                             { coordsName, coordinates },
                         },
        {} };
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelMetadata metadata;
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    grid.dumpModelState(state, metadata, filename, true);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    lun = std::fopen(filename.c_str(), "r");
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    REQUIRE(lun != NULL);
    REQUIRE(std::fclose(lun) == 0);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    // Reset the array dimensions to make sure that the read function gets them correct
    ModelArray::setDimension(ModelArray::Dimension::X, 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::Y, 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::Z, 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::XCG, 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ModelArray::setDimension(ModelArray::Dimension::YCG, 1);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    // In the full model numbers of DG components are set at compile time, so they are not reset
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    REQUIRE(ModelArray::nComponents(ModelArray::Type::VERTEX) == ModelArray::nCoords);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ParametricGrid gridIn;
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    ParaGridIO* readIO = new ParaGridIO(gridIn);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    gridIn.setIO(readIO);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelState ms = gridIn.getModelState(filename);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[0] == nx);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[1] == ny);
    REQUIRE(ModelArray::dimensions(ModelArray::Type::Z)[2] == nz);

    REQUIRE(ms.data.size() == state.data.size());
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelArray& ticeRef = ms.data.at(ticeName);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    REQUIRE(ModelArray::nDimensions(ModelArray::Type::Z) == 3);
    REQUIRE(ticeRef.getType() == ModelArray::Type::Z);
    REQUIRE(ticeRef.nDimensions() == 3);
    REQUIRE(ticeRef.dimensions()[0] == nx);
    REQUIRE(ticeRef.dimensions()[1] == ny);
    REQUIRE(ticeRef.dimensions()[2] == nz);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelArray& hiceRef = ms.data.at(hiceName);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    REQUIRE(hiceRef.nDimensions() == 2);
    REQUIRE(hiceRef.dimensions()[0] == nx);
    REQUIRE(hiceRef.dimensions()[1] == ny);
    REQUIRE(ModelArray::nComponents(ModelArray::Type::DG) == DG);
    REQUIRE(hiceRef.nComponents() == DG);
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    REQUIRE(ticeRef(12, 14, 1) == tice(12, 14, 1));
    std::cerr << "PGt checkpoint " << ++c << std::endl;

    ModelArray& coordRef = ms.data.at(coordsName);
    std::cerr << "PGt checkpoint " << ++c << std::endl;
    REQUIRE(coordRef.nDimensions() == 2);
    REQUIRE(coordRef.nComponents() == 2);
    REQUIRE(coordRef.components({ 12, 13 })[0] - coordRef.components({ 11, 13 })[0] == scale);
    REQUIRE(coordRef.components({ 12, 13 })[1] - coordRef.components({ 12, 12 })[1] == scale);

    std::remove(filename.c_str());
    std::cerr << "PG test 1 done" << std::endl;
}

TEST_CASE("Write a diagnostic ParaGrid file", "[ParaGridIO]")
{
    Module::setImplementation<IStructure>("ParametricGrid");

    REQUIRE(Module::getImplementation<IStructure>().structureType() == "parametric_rectangular");

    std::FILE* lun = std::fopen(diagFile.c_str(), "r");
    if (lun != NULL) {
        std::remove(diagFile.c_str());
    }
    std::fclose(lun);

    ParametricGrid grid;
    ParaGridIO* pio = new ParaGridIO(grid);
    grid.setIO(pio);

    // Set the dimension lengths
    size_t nx = 30;
    size_t ny = 20;
    size_t nz = 3;
    size_t nxcg = CG * nx + 1;
    size_t nycg = CG * ny + 1;

    double yFactor = 0.01;
    double xFactor = 0.0001;

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::Z, nz);
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
//                = (i - nx / 2) * (i - nx / 2) + (j - ny / 2) * (j - ny / 2) > (nx * ny) ? 0 : 1;
            for (size_t d = 0; d < DG; ++d) {
                fractionalDG.components({ i, j })[d] = fractional(i, j) + d;
            }
        }
    }
    double prec = 1e-9;
    REQUIRE(fractional(12, 12) - fractional(11, 12) == Approx(xFactor).epsilon(prec));
    REQUIRE(fractional(12, 12) - fractional(12, 11) == Approx(yFactor).epsilon(prec));

    REQUIRE(fractionalDG(12, 12) - fractionalDG(11, 12) == Approx(xFactor).epsilon(prec));
    REQUIRE(fractionalDG(12, 12) - fractionalDG(12, 11) == Approx(yFactor).epsilon(prec));


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

    ModelState state = { {
                             { maskName, mask },
                             { hiceName, hice },
                             { ciceName, cice },
                             { coordsName, coordinates },
                         },
        {} };

    ModelMetadata metadata;
    metadata.setTime(TimePoint("2000-01-01T00:00:00Z"));

    grid.dumpModelState(state, metadata, diagFile, false);

    for (int t = 1; t < 5; ++t) {
        hice += 100;
        cice += 100;
        state = { {
                      { hiceName, hice },
                      { ciceName, cice },
                  },
            {} };

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
    // TODO implement ParametricGrid in the IStructureModule
    REQUIRE(structureType == grid.structureType());

    // TODO test metadata

    // test data
    REQUIRE(dataGrp.getVarCount() == 4);
    netCDF::NcVar hiceVar = dataGrp.getVar(hiceName);
    netCDF::NcVar ciceVar = dataGrp.getVar(ciceName);
    netCDF::NcVar coordVar = dataGrp.getVar(coordsName);
    netCDF::NcVar maskVar = dataGrp.getVar(maskName);

    // hice
    REQUIRE(hiceVar.getDimCount() == 4);


    ncFile.close();

    std::remove(diagFile.c_str());

}

}
