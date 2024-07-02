/*!
 * @file    XiosDomain_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    27 June 2024
 * @brief   Tests for XIOS domains
 * @details
 * This test is designed to test domain functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

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
    REQUIRE_FALSE(xios_handler.isDefinedDomainType(domainId));
    const std::string domainType = { "rectilinear" };
    xios_handler.setDomainType(domainId, domainType);
    REQUIRE(xios_handler.isDefinedDomainType(domainId));
    REQUIRE(xios_handler.getDomainType(domainId) == domainType);
    // Global longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    const size_t ni_glo = 60;
    xios_handler.setDomainGlobalLongitudeSize(domainId, ni_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLongitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLongitudeSize(domainId) == ni_glo);
    // Global latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    const size_t nj_glo = 20;
    xios_handler.setDomainGlobalLatitudeSize(domainId, nj_glo);
    REQUIRE(xios_handler.isDefinedDomainGlobalLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainGlobalLatitudeSize(domainId) == nj_glo);
    // Local longitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeSize(domainId));
    const size_t ni = ni_glo / size;
    xios_handler.setDomainLongitudeSize(domainId, ni);
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLongitudeSize(domainId) == ni);
    // Local latitude size
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    const size_t nj = nj_glo;
    xios_handler.setDomainLatitudeSize(domainId, nj);
    REQUIRE(xios_handler.isDefinedDomainLatitudeSize(domainId));
    REQUIRE(xios_handler.getDomainLatitudeSize(domainId) == nj);
    // Local longitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    const size_t startLon = ni * rank;
    xios_handler.setDomainLongitudeStart(domainId, startLon);
    REQUIRE(xios_handler.isDefinedDomainLongitudeStart(domainId));
    REQUIRE(xios_handler.getDomainLongitudeStart(domainId) == startLon);
    // Local latitude start
    REQUIRE_FALSE(xios_handler.isDefinedDomainLatitudeStart(domainId));
    const size_t startLat = 0;
    xios_handler.setDomainLatitudeStart(domainId, startLat);
    REQUIRE(xios_handler.isDefinedDomainLatitudeStart(domainId));
    REQUIRE(xios_handler.getDomainLatitudeStart(domainId) == startLat);
    // Local longitude values
    REQUIRE_FALSE(xios_handler.areDefinedDomainLongitudeValues(domainId));
    std::vector<double> vecLon(ni);
    for (size_t i = 0; i < ni; i++) {
        vecLon[i] = -180 + (rank * ni * i) * 360 / ni_glo;
    }
    xios_handler.setDomainLongitudeValues(domainId, vecLon);
    REQUIRE(xios_handler.areDefinedDomainLongitudeValues(domainId));
    std::vector<double> vecLonOut = xios_handler.getDomainLongitudeValues(domainId);
    for (size_t i = 0; i < ni; i++) {
        REQUIRE(vecLonOut[i] == doctest::Approx(vecLon[i]));
    }
    // Local latitude values
    REQUIRE_FALSE(xios_handler.areDefinedDomainLatitudeValues(domainId));
    std::vector<double> vecLat(nj);
    for (size_t j = 0; j < nj; j++) {
        vecLat[j] = -90 + j * 180 / nj_glo;
    }
    xios_handler.setDomainLatitudeValues(domainId, vecLat);
    REQUIRE(xios_handler.areDefinedDomainLatitudeValues(domainId));
    std::vector<double> vecLatOut = xios_handler.getDomainLatitudeValues(domainId);
    for (size_t j = 0; j < nj; j++) {
        REQUIRE(vecLatOut[j] == doctest::Approx(vecLat[j]));
    }

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
