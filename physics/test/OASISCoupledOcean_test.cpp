/*!
 * @file OASISCoupledOcean_test.cpp
 *
 * @date 15 Feb 2025
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

    double sstIn = -1.84;
    double sssIn = 28.0;
    // double mldIn = 14.8;
    double uIn = -0.14;
    double vIn = 0.71;
    double sshIn = 14.8;

    HField cice(ModelArray::Type::H);
    cice = 0.8;
    ModelComponent::getStore().registerArray(Shared::C_ICE_DG, &cice, RO);

    HField tauXIO(ModelArray::Type::H);
    tauXIO = -3e-2;
    ModelComponent::getStore().registerArray(Protected::IO_STRESS_X, &tauXIO);

    HField tauYIO(ModelArray::Type::H);
    tauYIO = 4e-2;
    ModelComponent::getStore().registerArray(Protected::IO_STRESS_Y, &tauYIO);

    HField newIce(ModelArray::Type::H);
    newIce = 4e-2;
    ModelComponent::getStore().registerArray(Shared::NEW_ICE, &newIce, RW);

    HField deltaHice(ModelArray::Type::H);
    deltaHice = 1e-4;
    ModelComponent::getStore().registerArray(Shared::DELTA_HICE, &deltaHice, RW);

    HField deltaSmelt(ModelArray::Type::H);
    deltaSmelt = 1e-4;
    ModelComponent::getStore().registerArray(Shared::HSNOW_MELT, &deltaSmelt, RW);

    HField qow(ModelArray::Type::H);
    qow = 100;
    ModelComponent::getStore().registerArray(Shared::Q_OW, &qow, RW);

    HField tauXOW(ModelArray::Type::H);
    tauXOW = tauXIO;
    ModelComponent::getStore().registerArray(Shared::OW_STRESS_X, &tauXOW, RW);

    HField tauYOW(ModelArray::Type::H);
    tauYOW = tauYIO;
    ModelComponent::getStore().registerArray(Shared::OW_STRESS_Y, &tauYOW, RW);

    OASISCoupledOcean ocpl;
    ModelMetadata metadata;
    const std::vector<int> partInfo = { OASIS_Serial, 1, 1 };
    OASIS_CHECK_ERR(oasis_c_def_partition(&metadata.OASISPartitionId, OASIS_Serial_Params,
        &partInfo[0], OASIS_No_Gsize, compName.c_str()));

    ocpl.setData(ModelState::DataMap());
    ocpl.configure();
    ocpl.setMetadata(metadata);
    OASIS_CHECK_ERR(oasis_c_enddef());

    TimePoint t1("2000-01-01T00:00:00Z");
    TimestepTime tst = { t1, Duration(600) };

    ocpl.updateBefore(tst);

    ModelArrayRef<Protected::SST> sst(ModelComponent::getStore());
    ModelArrayRef<Protected::SSS> sss(ModelComponent::getStore());
    ModelArrayRef<Protected::OCEAN_U> u(ModelComponent::getStore());
    ModelArrayRef<Protected::OCEAN_V> v(ModelComponent::getStore());
    ModelArrayRef<Protected::SSH> ssh(ModelComponent::getStore());
    /*
    std::cout << "Received SST at time " << ocpl.OASISTime << ":  " << sst[0] << std::endl ;
    std::cout << "Received SSS at time " << ocpl.OASISTime << ":  " << sss[0] << std::endl ;
    std::cout << "Received OCEAN_U at time " << ocpl.OASISTime << ":  " << u[0] << std::endl ;
    std::cout << "Received OCEAN_V at time " << ocpl.OASISTime << ":  " << v[0] << std::endl ;
    std::cout << "Received OCEAN_V at time " << ocpl.OASISTime << ":  " << ssh[0] << std::endl ;
    */
    REQUIRE(sst[0] == sstIn);
    REQUIRE(sss[0] == sssIn);
    REQUIRE(u[0] == uIn);
    REQUIRE(v[0] == vIn);
    REQUIRE(ssh[0] == sshIn);
    // REQUIRE(mld[0] == mldIn);

    ModelArrayRef<Shared::Q_SW_BASE, RW> qSwBase(ModelComponent::getStore());
    qSwBase[0] = -2.;

    ModelArrayRef<Shared::Q_SW_OW, RW> qswow(ModelComponent::getStore());
    qswow[0] = -10.;

    ocpl.updateAfter(tst);

    /* The OASIS output file should contain:
     * I_taux: -0.03
     * I_tauy: 0.04
     * I_taumod: 0.05
     * I_fwflux: 0.0609935
     * I_rsso: 3.6
     * I_rsnos: 1651.14
     * I_sfi: 0.000306278
     * I_conc: 0.8
     */

    OASIS_CHECK_ERR(oasis_c_terminate());
}
TEST_SUITE_END();
}
