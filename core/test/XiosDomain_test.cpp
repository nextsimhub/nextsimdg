/*!
 * @file    XiosDomain_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    26 July 2024
 * @brief   Tests for XIOS domains
 * @details
 * This test is designed to test domain functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Configurator.hpp"
#include "include/Xios.hpp"

#include <iostream>

namespace Nextsim {

/*!
 * TestXiosDomin
 *
 * This function tests the domain functionality of the C++ interface for XIOS. It
 * needs to be run with 2 ranks i.e.,
 *
 * `mpirun -n 2 ./testXiosDomain_MPI2`
 *
 */
MPI_TEST_CASE("TestXiosDomain", 2)
{

    // Enable XIOS in the 'config'
    Configurator::clearStreams();
    std::stringstream config;
    config << "[xios]" << std::endl << "enable = true" << std::endl;
    std::unique_ptr<std::istream> pcstream(new std::stringstream(config.str()));
    Configurator::addStream(std::move(pcstream));

    // Initialize an Xios instance called xios_handler
    Xios xios_handler;
    REQUIRE(xios_handler.isInitialized());
    const size_t size = xios_handler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xios_handler.getClientMPIRank();

    // Set timestep as a minimum
    xios_handler.setCalendarTimestep(Duration("P0-0T01:00:00"));

    // --- Tests for domain API
    const std::string domainId = { "domain_A" };
    xios_handler.createDomain(domainId);
    // Domain type
    const std::string domainType = { "rectilinear" };
    xios_handler.setDomainType(domainId, domainType);
    REQUIRE(xios_handler.getDomainType(domainId) == domainType);
    // Global longitude size
    const size_t ni_glo = 60;
    xios_handler.setDomainGlobalXSize(domainId, ni_glo);
    REQUIRE(xios_handler.getDomainGlobalXSize(domainId) == ni_glo);
    // Global latitude size
    const size_t nj_glo = 20;
    xios_handler.setDomainGlobalYSize(domainId, nj_glo);
    REQUIRE(xios_handler.getDomainGlobalYSize(domainId) == nj_glo);
    // Local longitude size
    const size_t ni = ni_glo / size;
    xios_handler.setDomainLocalXSize(domainId, ni);
    REQUIRE(xios_handler.getDomainLocalXSize(domainId) == ni);
    // Local latitude size
    const size_t nj = nj_glo;
    xios_handler.setDomainLocalYSize(domainId, nj);
    REQUIRE(xios_handler.getDomainLocalYSize(domainId) == nj);
    // Local longitude start
    const size_t startLon = ni * rank;
    xios_handler.setDomainLocalXStart(domainId, startLon);
    REQUIRE(xios_handler.getDomainLocalXStart(domainId) == startLon);
    // Local latitude start
    const size_t startLat = 0;
    xios_handler.setDomainLocalYStart(domainId, startLat);
    REQUIRE(xios_handler.getDomainLocalYStart(domainId) == startLat);
    // Local longitude values
    std::vector<double> vx(ni);
    for (size_t i {}; i < ni; i++) {
        vx[i] = -180 + (rank * ni * i) * 360 / ni_glo;
    }
    xios_handler.setDomainLocalXValues(domainId, vx);
    std::vector<double> vxOut = xios_handler.getDomainLocalXValues(domainId);
    for (size_t i {}; i < ni; i++) {
        REQUIRE(vxOut[i] == doctest::Approx(vx[i]));
    }
    // Local latitude values
    std::vector<double> vy(nj);
    for (size_t j {}; j < nj; j++) {
        vy[j] = -90 + j * 180 / nj_glo;
    }
    xios_handler.setDomainLocalYValues(domainId, vy);
    std::vector<double> vyOut = xios_handler.getDomainLocalYValues(domainId);
    for (size_t j {}; j < nj; j++) {
        REQUIRE(vyOut[j] == doctest::Approx(vy[j]));
    }

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
