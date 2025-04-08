/*!
 * @file    XiosDomain_test.cpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    08 Apr 2025
 * @brief   Tests for XIOS domains
 * @details
 * This test is designed to test domain functionality of the C++ interface
 * for XIOS.
 *
 */
#include <doctest/extensions/doctest_mpi.h>
#undef INFO

#include "StructureModule/include/ParametricGrid.hpp"
#include "include/Xios.hpp"

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
    enableXios();

    // Get the Xios singleton instance and check it's initialized
    Xios& xiosHandler = Xios::getInstance();
    REQUIRE(xiosHandler.isInitialized());
    const size_t size = xiosHandler.getClientMPISize();
    REQUIRE(size == 2);
    const size_t rank = xiosHandler.getClientMPIRank();

    // --- Tests for domain API
    const std::string domainId = "domain_A";
    REQUIRE_THROWS_WITH(xiosHandler.getDomainType(domainId), "Xios: Undefined domain 'domain_A'");
    xiosHandler.createDomain(domainId);
    REQUIRE_THROWS_WITH(
        xiosHandler.createDomain(domainId), "Xios: Domain 'domain_A' already exists");
    // Domain type
    REQUIRE_THROWS_WITH(
        xiosHandler.getDomainType(domainId), "Xios: Undefined type for domain 'domain_A'");
    const std::string domainType = "rectilinear";
    xiosHandler.setDomainType(domainId, domainType);
    REQUIRE(xiosHandler.getDomainType(domainId) == domainType);
    // Global number of points in x-direction
    REQUIRE_THROWS_WITH(xiosHandler.getDomainGlobalXSize(domainId),
        "Xios: Undefined global x-size for domain 'domain_A'");
    const size_t nx_glo = 4;
    xiosHandler.setDomainGlobalXSize(domainId, nx_glo);
    REQUIRE(xiosHandler.getDomainGlobalXSize(domainId) == nx_glo);
    // Global number of points in y-direction
    REQUIRE_THROWS_WITH(xiosHandler.getDomainGlobalYSize(domainId),
        "Xios: Undefined global y-size for domain 'domain_A'");
    const size_t ny_glo = 2;
    xiosHandler.setDomainGlobalYSize(domainId, ny_glo);
    REQUIRE(xiosHandler.getDomainGlobalYSize(domainId) == ny_glo);
    // Local number of points in x-direction
    REQUIRE_THROWS_WITH(xiosHandler.getDomainLocalXSize(domainId),
        "Xios: Undefined local x-size for domain 'domain_A'");
    const size_t nx = nx_glo / size;
    xiosHandler.setDomainLocalXSize(domainId, nx);
    REQUIRE(xiosHandler.getDomainLocalXSize(domainId) == nx);
    // Local number of points in y-direction
    REQUIRE_THROWS_WITH(xiosHandler.getDomainLocalYSize(domainId),
        "Xios: Undefined local y-size for domain 'domain_A'");
    const size_t ny = ny_glo;
    xiosHandler.setDomainLocalYSize(domainId, ny);
    REQUIRE(xiosHandler.getDomainLocalYSize(domainId) == ny);
    // Local starting x-index
    REQUIRE_THROWS_WITH(xiosHandler.getDomainLocalXStart(domainId),
        "Xios: Undefined local starting x-index for domain 'domain_A'");
    const size_t x0 = nx * rank;
    xiosHandler.setDomainLocalXStart(domainId, x0);
    REQUIRE(xiosHandler.getDomainLocalXStart(domainId) == x0);
    // Local starting y-index
    REQUIRE_THROWS_WITH(xiosHandler.getDomainLocalYStart(domainId),
        "Xios: Undefined local starting y-index for domain 'domain_A'");
    const size_t y0 = 0;
    xiosHandler.setDomainLocalYStart(domainId, y0);
    REQUIRE(xiosHandler.getDomainLocalYStart(domainId) == y0);
    // Local x-values
    REQUIRE_THROWS_WITH(xiosHandler.getDomainLocalXValues(domainId),
        "Xios: Undefined local x-values for domain 'domain_A'");
    REQUIRE_THROWS_AS(xiosHandler.getDomainLocalXValues(domainId), std::runtime_error);
    std::vector<double> vx { -1.0 + rank, -0.5 + rank };
    xiosHandler.setDomainLocalXValues(domainId, vx);
    std::vector<double> vxOut = xiosHandler.getDomainLocalXValues(domainId);
    for (size_t i = 0; i < nx; i++) {
        REQUIRE(vxOut[i] == doctest::Approx(vx[i]));
    }
    // Local y-values
    REQUIRE_THROWS_WITH(xiosHandler.getDomainLocalYValues(domainId),
        "Xios: Undefined local y-values for domain 'domain_A'");
    std::vector<double> vy { -1.0, 1.0 };
    xiosHandler.setDomainLocalYValues(domainId, vy);
    std::vector<double> vyOut = xiosHandler.getDomainLocalYValues(domainId);
    for (size_t j = 0; j < ny; j++) {
        REQUIRE(vyOut[j] == doctest::Approx(vy[j]));
    }

    xiosHandler.close_context_definition();
    xiosHandler.context_finalize();
}
}
