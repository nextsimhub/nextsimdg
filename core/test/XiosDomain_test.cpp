/*!
 * @file    XiosDomain_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    21 August 2024
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
 * TestXiosDomain
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
    // Global number of points in x-direction
    const size_t nx_glo = 4;
    xios_handler.setDomainGlobalXSize(domainId, nx_glo);
    REQUIRE(xios_handler.getDomainGlobalXSize(domainId) == nx_glo);
    // Global number of points in y-direction
    const size_t ny_glo = 2;
    xios_handler.setDomainGlobalYSize(domainId, ny_glo);
    REQUIRE(xios_handler.getDomainGlobalYSize(domainId) == ny_glo);
    // Local number of points in x-direction
    const size_t nx = nx_glo / size;
    xios_handler.setDomainLocalXSize(domainId, nx);
    REQUIRE(xios_handler.getDomainLocalXSize(domainId) == nx);
    // Local number of points in y-direction
    const size_t ny = ny_glo;
    xios_handler.setDomainLocalYSize(domainId, ny);
    REQUIRE(xios_handler.getDomainLocalYSize(domainId) == ny);
    // Local starting x-index
    const size_t x0 = nx * rank;
    xios_handler.setDomainLocalXStart(domainId, x0);
    REQUIRE(xios_handler.getDomainLocalXStart(domainId) == x0);
    // Local starting y-index
    const size_t y0 = 0;
    xios_handler.setDomainLocalYStart(domainId, y0);
    REQUIRE(xios_handler.getDomainLocalYStart(domainId) == y0);
    // Local x-values
    std::vector<double> vx { -1.0 + rank, -0.5 + rank };
    xios_handler.setDomainLocalXValues(domainId, vx);
    std::vector<double> vxOut = xios_handler.getDomainLocalXValues(domainId);
    for (size_t i = 0; i < nx; i++) {
        REQUIRE(vxOut[i] == doctest::Approx(vx[i]));
    }
    // Local y-values
    std::vector<double> vy { -1.0, 1.0 };
    xios_handler.setDomainLocalYValues(domainId, vy);
    std::vector<double> vyOut = xios_handler.getDomainLocalYValues(domainId);
    for (size_t j = 0; j < ny; j++) {
        REQUIRE(vyOut[j] == doctest::Approx(vy[j]));
    }

    xios_handler.close_context_definition();
    xios_handler.context_finalize();
}
}
