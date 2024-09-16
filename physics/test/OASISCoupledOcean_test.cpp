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
MPI_TEST_CASE("OASIS init put and get", 1)
{
    MPI_Comm modelCommunicator;
    int compID; // Not actually used. Only useful for debugging
    const std::string compName = "nextsim"; // Not useful for any setups we have in mind
    OASIS_CHECK_ERR(oasis_c_init_comp(&compID, compName.c_str(), OASIS_COUPLED));
    OASIS_CHECK_ERR(oasis_c_get_localcomm(&modelCommunicator));

    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModelArray::setDimensions(ModelArray::Type::Z, { 1, 1, 1 });

    double sstIn = -1.84;
    double sssIn = 28.0;
    //double mldIn = 14.8;
    double uIn = -0.14;
    double vIn = 0.71;

    HField cice(ModelArray::Type::H);
    cice = 1.0;
    ModelComponent::getStore().registerArray(Protected::C_ICE, &cice, RO);
    OASISCoupledOcean ocpl;
    ModelMetadata metadata;
    const std::vector<int> partInfo = { OASIS_Serial, 1, 1 };
    OASIS_CHECK_ERR(oasis_c_def_partition(
        &metadata.OASISPartitionId, OASIS_Serial_Params, &partInfo[0], OASIS_No_Gsize, compName.c_str()));

    ocpl.setData(ModelState::DataMap());
    ocpl.configure();
    ocpl.setMetadata(metadata);
    OASIS_CHECK_ERR(oasis_c_enddef());

    ocpl.updateBefore(TimestepTime());
    ModelArrayRef<Protected::SST> sst(ModelComponent::getStore());
    ModelArrayRef<Protected::SSS> sss(ModelComponent::getStore());
    ModelArrayRef<Protected::OCEAN_U> u(ModelComponent::getStore());
    ModelArrayRef<Protected::OCEAN_V> v(ModelComponent::getStore());
    std::cout << "Received SST at time " << ocpl.OASISTime << ":  " << sst[0] << std::endl ;
    std::cout << "Received SSS at time " << ocpl.OASISTime << ":  " << sss[0] << std::endl ;
    std::cout << "Received OCEAN_U at time " << ocpl.OASISTime << ":  " << u[0] << std::endl ;
    std::cout << "Received OCEAN_V at time " << ocpl.OASISTime << ":  " << v[0] << std::endl ;
    REQUIRE(sst[0] == sstIn);
    REQUIRE(sss[0] == sssIn);
    REQUIRE(u[0] == uIn);
    REQUIRE(v[0] == vIn);
    //REQUIRE(mld[0] == mldIn);

    ocpl.updateAfter(TimestepTime());

    OASIS_CHECK_ERR(oasis_c_terminate());
}
TEST_SUITE_END();
}
