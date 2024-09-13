/*!
 * @file OASISCoupledOcean_test.cpp
 *
 * @date 13 Sep 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#include <doctest/extensions/doctest_mpi.h>

#include "include/OASISCoupledOcean.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("OASISCoupledOcean");
MPI_TEST_CASE("OASIS init put and get", 2)
{
    MPI_Comm modelCommunicator;
    int compID; // Not actually used. Only useful for debugging
    const std::string compName = "nextsim"; // Not useful for any setups we have in mind
    OASIS_CHECK_ERR(oasis_c_init_comp(&compID, compName.c_str(), OASIS_COUPLED));
    OASIS_CHECK_ERR(oasis_c_get_localcomm(&modelCommunicator));

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    HField cice(ModelArray::Type::H);
    cice = 1.0;
    ModelComponent::getStore().registerArray(Protected::C_ICE, &cice, RO);
    OASISCoupledOcean ocpl;
    ModelMetadata metadata;

    ocpl.setData(ModelState::DataMap());
    ocpl.setMetadata(metadata);
    ocpl.updateBefore(TimestepTime());
    ocpl.updateAfter(TimestepTime());

    OASIS_CHECK_ERR(oasis_c_terminate());
}
TEST_SUITE_END();
}