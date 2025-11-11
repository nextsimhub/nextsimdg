/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ConstantOceanBoundary.hpp"

#include "include/ModelState.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("ConstantOceanBoundary");
TEST_CASE("ConstantOcean Qio calculation")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor(
        ModelComponent::getStore(), RO, ModelArray::Type::H);
    HField& cice = ciceAccessor.getHostRW();
    cice = 1.0; // Need some ice if Qio is to be calculated
    ConstantOceanBoundary cob;

    cob.setData(ModelState::DataMap());
    cob.updateBefore(TimestepTime());
    
    ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor(ModelComponent::getStore());
    const HField& qio = qioAccessor.getHostRO();

    REQUIRE(qio[0] != 0.);
}
TEST_SUITE_END();
}
