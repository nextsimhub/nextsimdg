/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief   Test halo exchange class for Closed Boundaries (CB)
 * @details
 *
 * Test halo exchange class for Closed Boundaries (CB) on following fields:
 * - HField
 * - VertexField
 * - DGField
 * - DGVector
 */

#include <doctest/extensions/doctest_mpi.h>

#include "ModelMPI.hpp"
#include "ModelMetadata.hpp"
#include "cgVector.hpp"
#include "include/DGModelArray.hpp"
#include "include/Halo.hpp"
#include "include/Interpolations.hpp"

namespace Nextsim {

const std::string testFilesDir = TEST_FILES_DIR;
const std::string file = testFilesDir + "/partition_metadata_3_cb.nc";

static const int DG = 3;
static const int CGDEGREE = 2;
static const bool debug = false;

// reference data for each process
static std::array<std::vector<double>, 3> refDataAllProcs = {
    std::vector<double>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 100, 110, 120, 130, 140, 150, 160, 370, 0,
        101, 111, 121, 131, 141, 151, 161, 371, 0, 102, 112, 122, 132, 142, 152, 162, 372, 0, 103,
        113, 123, 133, 143, 153, 163, 373, 0, 204, 214, 224, 234, 244, 254, 264, 0 }),
    std::vector<double>({ 0, 103, 113, 123, 133, 143, 153, 163, 0, 0, 204, 214, 224, 234, 244, 254,
        264, 374, 0, 205, 215, 225, 235, 245, 255, 265, 375, 0, 206, 216, 226, 236, 246, 256, 266,
        376, 0, 207, 217, 227, 237, 247, 257, 267, 377, 0, 208, 218, 228, 238, 248, 258, 268, 378,
        0, 0, 0, 0, 0, 0, 0, 0, 0 }),
    std::vector<double>({ 0, 0, 0, 0, 0, 160, 370, 380, 390, 0, 161, 371, 381, 391, 0, 162, 372,
        382, 392, 0, 163, 373, 383, 393, 0, 264, 374, 384, 394, 0, 265, 375, 385, 395, 0, 266, 376,
        386, 396, 0, 267, 377, 387, 397, 0, 268, 378, 388, 398, 0, 0, 0, 0, 0, 0 }),
};

void initializeTestData(ModelArray::DataType& data, size_t localNx, size_t localNy, size_t offsetX,
    size_t offsetY, size_t numComps, int test_rank)
{
    // initialize with test data
    for (size_t j = 0; j < localNy; ++j) {
        for (size_t i = 0; i < localNx; ++i) {
            for (size_t d = 0; d < numComps; ++d) {
                data(d + i * numComps + j * localNx * numComps)
                    = (d + 1) * 1000 + (test_rank + 1) * 100 + (i + offsetX) * 10 + (j + offsetY);
            }
        }
    }
}

TEST_SUITE_BEGIN("Halo exchange tests");
MPI_TEST_CASE("test halo exchange on 3 proc grid", 3)
{
    // test for a 3 proc grid
    // ┌─┬─┐
    // │1│ │
    // ├─┤2│
    // │0│ │
    // └─┴─┘ (proc id)

    auto& modelMPI = ModelMPI::getInstance(test_comm);
    auto& metadata = ModelMetadata::getInstance(file);

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1, localNx + 1, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1, localNy + 1, offsetY);

    // create example 2D field on each process
    auto testData = ModelArray::HField();
    testData.resize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), testData.nComponents());
    initializeTestData(
        innerData, localNx - 2, localNy - 2, offsetX, offsetY, testData.nComponents(), test_rank);

    // populate inner block of modelarray
    halo.setInnerBlock(innerData, testData.getDataRef());
    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    // initialize reference data for each proc
    ModelArray::DataType refData;
    refData.resize(localNx * localNy, 1);
    refData
        = Eigen::Map<ModelArray::DataType>(refDataAllProcs[test_rank].data(), refData.size(), 1);

    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            if (refData(i + j * localNx) > 0) {
                REQUIRE(testData(i, j) == 1000. + refData(i + j * localNx));
            } else {
                REQUIRE(testData(i, j) == 0.);
            }
        }
    }

    VertexField coordinates = ModelArray::VertexField();
    coordinates.resize();
    coordinates = 0.;

    Halo haloVertex(coordinates);

    // Vetex coordinates
    auto localNxVertex = localNx + 1;
    auto localNyVertex = localNy + 1;
    for (size_t i = 0; i < localNxVertex - 2 * Halo::haloWidth; ++i) {
        for (size_t j = 0; j < localNyVertex - 2 * Halo::haloWidth; ++j) {
            double x = (i + offsetX) - 0.5 - float(nx) / 2;
            double y = (j + offsetY) - 0.5 - float(ny) / 2;
            coordinates.components({ i + Halo::haloWidth, j + Halo::haloWidth })[0] = x * 100.;
            coordinates.components({ i + Halo::haloWidth, j + Halo::haloWidth })[1] = y * 100.;
        }
    }

    haloVertex.exchangeHalos(coordinates.getDataRef());

    std::vector<double> refDataVertex;

    if (test_rank == 0) {
        refDataVertex = std::vector<double>({ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., -550., -500., -450., -500., -350., -500., -250.,
            -500., -150., -500., -50., -500., 50., -500., 150., -500., 250., -500., 0., 0., -550.,
            -400., -450., -400., -350., -400., -250., -400., -150., -400., -50., -400., 50., -400.,
            150., -400., 250., -400., 0., 0., -550., -300., -450., -300., -350., -300., -250.,
            -300., -150., -300., -50., -300., 50., -300., 150., -300., 250., -300., 0., 0., -550.,
            -200., -450., -200., -350., -200., -250., -200., -150., -200., -50., -200., 50., -200.,
            150., -200., 250., -200., 0., 0., -550., -100., -450., -100., -350., -100., -250.,
            -100., -150., -100., -50., -100., 50., -100., 150., -100., 250., -100., 0., 0., -550.,
            0., -450., 0., -350., 0., -250., 0., -150., 0., -50., 0., 50., 0., 150., 0., 0., 0. });
    }
    if (test_rank == 1) {
        refDataVertex = std::vector<double>({ 0., 0., -550., -200., -450., -200., -350., -200.,
            -250., -200., -150., -200., -50., -200., 50., -200., 150., -200., 0., 0., 0., 0., -550.,
            -100., -450., -100., -350., -100., -250., -100., -150., -100., -50., -100., 50., -100.,
            150., -100., 250., -100., 0., 0., -550., 0., -450., 0., -350., 0., -250., 0., -150., 0.,
            -50., 0., 50., 0., 150., 0., 250., 0., 0., 0., -550., 100., -450., 100., -350., 100.,
            -250., 100., -150., 100., -50., 100., 50., 100., 150., 100., 250., 100., 0., 0., -550.,
            200., -450., 200., -350., 200., -250., 200., -150., 200., -50., 200., 50., 200., 150.,
            200., 250., 200., 0., 0., -550., 300., -450., 300., -350., 300., -250., 300., -150.,
            300., -50., 300., 50., 300., 150., 300., 250., 300., 0., 0., -550., 400., -450., 400.,
            -350., 400., -250., 400., -150., 400., -50., 400., 50., 400., 150., 400., 250., 400.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. });
    }
    if (test_rank == 2) {
        refDataVertex = std::vector<double>({ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 50.,
            -500., 150., -500., 250., -500., 350., -500., 450., -500., 0., 0., 50., -400., 150.,
            -400., 250., -400., 350., -400., 450., -400., 0., 0., 50., -300., 150., -300., 250.,
            -300., 350., -300., 450., -300., 0., 0., 50., -200., 150., -200., 250., -200., 350.,
            -200., 450., -200., 0., 0., 50., -100., 150., -100., 250., -100., 350., -100., 450.,
            -100., 0., 0., 50., 0., 150., 0., 250., 0., 350., 0., 450., 0., 0., 0., 50., 100., 150.,
            100., 250., 100., 350., 100., 450., 100., 0., 0., 50., 200., 150., 200., 250., 200.,
            350., 200., 450., 200., 0., 0., 50., 300., 150., 300., 250., 300., 350., 300., 450.,
            300., 0., 0., 50., 400., 150., 400., 250., 400., 350., 400., 450., 400., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0. });
    }

    auto coordRef = coordinates.getDataRef();
    for (size_t j = 0; j < localNyVertex; j++) {
        for (size_t i = 0; i < localNxVertex; i++) {
            for (size_t coord = 0; coord < 2; coord++) {
                REQUIRE(coordRef(i + j * localNxVertex, coord)
                    == refDataVertex[coord + i * 2 + j * localNxVertex * 2]);
            }
        }
    }
}

MPI_TEST_CASE("DGField", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    // create example 2D field on each process
    auto testData = ModelArray::DGField();
    testData.resize();
    testData = 0.;

    // create halo for testData model array
    Halo halo(testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), testData.nComponents());
    initializeTestData(
        innerData, localNx - 2, localNy - 2, offsetX, offsetY, testData.nComponents(), test_rank);

    // populate inner block of modelarray
    halo.setInnerBlock(innerData, testData.getDataRef());
    // exchange halos
    halo.exchangeHalos(testData.getDataRef());

    // initialize reference data for each proc
    ModelArray::DataType refData;
    refData.resize(localNx * localNy, 1);
    refData
        = Eigen::Map<ModelArray::DataType>(refDataAllProcs[test_rank].data(), refData.size(), 1);

    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            for (size_t d = 0; d < DG; ++d) {
                if (refData(i + j * localNx) > 0) {
                    REQUIRE(testData.components({ i, j })[d]
                        == (d + 1) * 1000. + refData(i + j * localNx));
                } else {
                    REQUIRE(testData.components({ i, j })[d] == 0.);
                }
            }
        }
    }
}

MPI_TEST_CASE("DGVector", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ParametricMesh smesh(CARTESIAN);
    smesh.nx = localNx;
    smesh.ny = localNy;
    smesh.nelements = localNx * localNy;
    smesh.nnodes = (localNx + 1) * (localNy + 1);
    smesh.vertices.resize(smesh.nnodes, Eigen::NoChange);
    for (size_t i = 0; i < localNx + 1; ++i) {
        for (size_t j = 0; j < localNy + 1; ++j) {
            smesh.vertices(i * (localNy + 1) + j, 0) = i;
            smesh.vertices(i * (localNy + 1) + j, 1) = j;
        }
    }

    // create example DGVector
    DGVector<DG> testData(smesh);
    testData.zero();

    // create halo for testData model array
    Halo halo(testData);

    // create and allocate temporary Eigen array
    ModelArray::DataType innerData;
    innerData.resize(halo.getInnerSize(), DG);
    initializeTestData(innerData, localNx - 2, localNy - 2, offsetX, offsetY, DG, test_rank);

    // populate inner block of modelarray
    halo.setInnerBlock(innerData, testData);

    // exchange halos
    halo.exchangeHalos(testData);

    // initialize reference data for each proc
    ModelArray::DataType refData;
    refData.resize(localNx * localNy, 1);
    refData
        = Eigen::Map<ModelArray::DataType>(refDataAllProcs[test_rank].data(), refData.size(), 1);

    for (size_t j = 0; j < localNy; j++) {
        for (size_t i = 0; i < localNx; i++) {
            for (size_t d = 0; d < DG; ++d) {
                if (refData(i + j * localNx) > 0) {
                    REQUIRE(
                        testData(i + j * localNx, d) == (d + 1) * 1000. + refData(i + j * localNx));
                } else {
                    REQUIRE(testData(i + j * localNx, d) == 0.);
                }
            }
        }
    }

    // initialize CGVector (after halo exchange on DGVector)
    CGVector<CGDEGREE> cgVector;
    cgVector.resize_by_mesh(smesh);
    Interpolations::DG2CG(smesh, cgVector, testData);

    // create halo for testData model array
    Halo haloCG(testData);
    // haloCG.exchangeHalos(cgVector);
}

MPI_TEST_CASE("CGVector", 3)
{
    auto& modelMPI = ModelMPI::getInstance();
    auto& metadata = ModelMetadata::getInstance();

    const size_t nx = metadata.getGlobalExtentX();
    const size_t ny = metadata.getGlobalExtentY();
    const size_t localNx = metadata.getLocalExtentX() + 2 * Halo::haloWidth;
    const size_t localNy = metadata.getLocalExtentY() + 2 * Halo::haloWidth;
    const size_t offsetX = metadata.getLocalCornerX();
    const size_t offsetY = metadata.getLocalCornerY();

    ModelArray::setDimension(ModelArray::Dimension::X, nx, localNx, offsetX);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny, localNy, offsetY);
    ModelArray::setNComponents(ModelArray::Type::DG, DG);

    ParametricMesh smesh(CARTESIAN);
    smesh.nx = localNx;
    smesh.ny = localNy;
    smesh.nelements = localNx * localNy;
    smesh.nnodes = (localNx + 1) * (localNy + 1);
    smesh.vertices.resize(smesh.nnodes, Eigen::NoChange);
    for (size_t i = 0; i < localNx + 1; ++i) {
        for (size_t j = 0; j < localNy + 1; ++j) {
            smesh.vertices(i * (localNy + 1) + j, 0) = i;
            smesh.vertices(i * (localNy + 1) + j, 1) = j;
        }
    }

    auto CGNx = localNx * CGDEGREE + 1;
    auto CGNy = localNy * CGDEGREE + 1;
    auto CGoffsetX = offsetX * CGDEGREE + 1;
    auto CGoffsetY = offsetY * CGDEGREE + 1;

    // create CGVector
    CGVector<CGDEGREE> cgVector;
    cgVector.resize_by_mesh(smesh);
    cgVector.zero();

    // Interpolations::DG2CG(smesh, cgVector, testData);

    // initialize data
    for (size_t j = CGDEGREE; j < CGNy - CGDEGREE; ++j) {
        for (size_t i = CGDEGREE; i < CGNx - CGDEGREE; ++i) {
            cgVector(i + j * CGNx)
                = ((i - CGDEGREE) + CGoffsetX) * 100 + ((j - CGDEGREE) + CGoffsetY);
        }
    }

    // create halo for testData model array
    Halo halo(cgVector);
    halo.exchangeHalos(cgVector);

    std::array<std::vector<double>, 3> cgRefDataAllProcs = {
        std::vector<double>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 101, 201, 301, 401, 501, 601, 701,
            801, 901, 1001, 1101, 1201, 1301, 1401, 1501, 1601, 1701, 0, 0, 102, 202, 302, 402, 502,
            602, 702, 802, 902, 1002, 1102, 1202, 1302, 1402, 1502, 1602, 1702, 0, 0, 103, 203, 303,
            403, 503, 603, 703, 803, 903, 1003, 1103, 1203, 1303, 1403, 1503, 1603, 1703, 0, 0, 104,
            204, 304, 404, 504, 604, 704, 804, 904, 1004, 1104, 1204, 1304, 1404, 1504, 1604, 1704,
            0, 0, 105, 205, 305, 405, 505, 605, 705, 805, 905, 1005, 1105, 1205, 1305, 1405, 1505,
            1605, 1705, 0, 0, 106, 206, 306, 406, 506, 606, 706, 806, 906, 1006, 1106, 1206, 1306,
            1406, 1506, 1606, 1706, 0, 0, 107, 207, 307, 407, 507, 607, 707, 807, 907, 1007, 1107,
            1207, 1307, 1407, 1507, 1607, 1707, 0, 0, 108, 208, 308, 408, 508, 608, 708, 808, 908,
            1008, 1108, 1208, 1308, 1408, 1508, 1608, 1708, 0, 0, 109, 209, 309, 409, 509, 609, 709,
            809, 909, 1009, 1109, 1209, 1309, 1409, 1509, 1609, 1709, 0, 0, 110, 210, 310, 410, 510,
            610, 710, 810, 910, 1010, 1110, 1210, 1310, 1410, 1510, 0, 0, 0, 0, 111, 211, 311, 411,
            511, 611, 711, 811, 911, 1011, 1111, 1211, 1311, 1411, 1511, 0, 0 }),
        std::vector<double>({ 0, 0, 107, 207, 307, 407, 507, 607, 707, 807, 907, 1007, 1107, 1207,
            1307, 1407, 1507, 0, 0, 0, 0, 108, 208, 308, 408, 508, 608, 708, 808, 908, 1008, 1108,
            1208, 1308, 1408, 1508, 0, 0, 0, 0, 109, 209, 309, 409, 509, 609, 709, 809, 909, 1009,
            1109, 1209, 1309, 1409, 1509, 1609, 1709, 0, 0, 110, 210, 310, 410, 510, 610, 710, 810,
            910, 1010, 1110, 1210, 1310, 1410, 1510, 1610, 1710, 0, 0, 111, 211, 311, 411, 511, 611,
            711, 811, 911, 1011, 1111, 1211, 1311, 1411, 1511, 1611, 1711, 0, 0, 112, 212, 312, 412,
            512, 612, 712, 812, 912, 1012, 1112, 1212, 1312, 1412, 1512, 1612, 1712, 0, 0, 113, 213,
            313, 413, 513, 613, 713, 813, 913, 1013, 1113, 1213, 1313, 1413, 1513, 1613, 1713, 0, 0,
            114, 214, 314, 414, 514, 614, 714, 814, 914, 1014, 1114, 1214, 1314, 1414, 1514, 1614,
            1714, 0, 0, 115, 215, 315, 415, 515, 615, 715, 815, 915, 1015, 1115, 1215, 1315, 1415,
            1515, 1615, 1715, 0, 0, 116, 216, 316, 416, 516, 616, 716, 816, 916, 1016, 1116, 1216,
            1316, 1416, 1516, 1616, 1716, 0, 0, 117, 217, 317, 417, 517, 617, 717, 817, 917, 1017,
            1117, 1217, 1317, 1417, 1517, 1617, 1717, 0, 0, 118, 218, 318, 418, 518, 618, 718, 818,
            918, 1018, 1118, 1218, 1318, 1418, 1518, 1618, 1718, 0, 0, 119, 219, 319, 419, 519, 619,
            719, 819, 919, 1019, 1119, 1219, 1319, 1419, 1519, 1619, 1719, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0 }),
        std::vector<double>({ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1301, 1401, 1501, 1601, 1701, 1801, 1901, 2001, 2101, 0, 0, 1302, 1402, 1502, 1602,
            1702, 1802, 1902, 2002, 2102, 0, 0, 1303, 1403, 1503, 1603, 1703, 1803, 1903, 2003,
            2103, 0, 0, 1304, 1404, 1504, 1604, 1704, 1804, 1904, 2004, 2104, 0, 0, 1305, 1405,
            1505, 1605, 1705, 1805, 1905, 2005, 2105, 0, 0, 1306, 1406, 1506, 1606, 1706, 1806,
            1906, 2006, 2106, 0, 0, 1307, 1407, 1507, 1607, 1707, 1807, 1907, 2007, 2107, 0, 0,
            1308, 1408, 1508, 1608, 1708, 1808, 1908, 2008, 2108, 0, 0, 1309, 1409, 1509, 1609,
            1709, 1809, 1909, 2009, 2109, 0, 0, 1310, 1410, 1510, 1610, 1710, 1810, 1910, 2010,
            2110, 0, 0, 1311, 1411, 1511, 1611, 1711, 1811, 1911, 2011, 2111, 0, 0, 1312, 1412,
            1512, 1612, 1712, 1812, 1912, 2012, 2112, 0, 0, 1313, 1413, 1513, 1613, 1713, 1813,
            1913, 2013, 2113, 0, 0, 1314, 1414, 1514, 1614, 1714, 1814, 1914, 2014, 2114, 0, 0,
            1315, 1415, 1515, 1615, 1715, 1815, 1915, 2015, 2115, 0, 0, 1316, 1416, 1516, 1616,
            1716, 1816, 1916, 2016, 2116, 0, 0, 1317, 1417, 1517, 1617, 1717, 1817, 1917, 2017,
            2117, 0, 0, 1318, 1418, 1518, 1618, 1718, 1818, 1918, 2018, 2118, 0, 0, 1319, 1419,
            1519, 1619, 1719, 1819, 1919, 2019, 2119, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0 }),
    };

    for (size_t j = 0; j < CGNy; j++) {
        for (size_t i = 0; i < CGNx; i++) {
            REQUIRE(cgVector(i + j * CGNx) == cgRefDataAllProcs[test_rank][i + j * CGNx]);
        }
    }
}
}
