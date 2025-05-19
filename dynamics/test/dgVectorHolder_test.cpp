/*!
 * @file DGVectorHolder_test.cpp
 *
 * @date Feb 4, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/dgVectorHolder.hpp"

namespace Nextsim {

COORDINATES CoordinateSystem = CARTESIAN;

TEST_SUITE_BEGIN("DGVectorHolder");
TEST_CASE("Cast to and from")
{
    static const int DG = 3;
    const size_t nx = 32;
    const size_t ny = 32;

    ParametricMesh smesh(CoordinateSystem);
    smesh.nx = nx;
    smesh.ny = ny;
    smesh.nelements = nx * ny;

    DGVector<DG>::EigenDGVector edgv(smesh.nelements, DG);
    double e0 = 3.;
    edgv.setConstant(e0);
    DGVectorHolder<DG> eHolder(edgv);
    REQUIRE(eHolder(0, 1) == e0);

    DGVector<DG> dgv;
    dgv.resize_by_mesh(smesh);
    double d0 = 4.;
    dgv.setConstant(d0);
    DGVectorHolder<DG> dHolder(dgv);
    REQUIRE(dHolder(0, 1) == d0);

    ModelArray::setDimension(ModelArray::Dimension::X, nx);
    ModelArray::setDimension(ModelArray::Dimension::Y, ny);
    ModelArray::setDimension(ModelArray::Dimension::DG, DG);

    ModelArray ma(ModelArray::Type::DG);
    ma.resize();
    double ma0 = 2.;
    ma = ma0;
    DGVectorHolder<DG> mHolder(ma);
    ModelArray::DataType& ema = static_cast<ModelArray::DataType&>(ma);
    REQUIRE(ma.components(0)[1] == ma0);
    REQUIRE(ema(0, 1) == ma0);
    REQUIRE(mHolder(0, 1) == ma0);

}
}

