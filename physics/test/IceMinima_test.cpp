/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/IceMinima.hpp"

namespace Nextsim {

class IColumnPhysics {
public:
    inline static void setMinima(double hMin, double cMin)
    {
        IceMinima::hMin = hMin;
        IceMinima::cMin = cMin;
    }
};

TEST_SUITE_BEGIN("IceMinima");
TEST_CASE("Set and retrieve ice minima values")
{
    // IceGrowth is a friend of IceMinima. Here is a class IceGrowth that
    // directly sets the minimum values without configuration.
    const double hMinTest = 0.1;
    const double cMinTest = 1e-6;

    IColumnPhysics::setMinima(hMinTest, cMinTest);
    REQUIRE(IceMinima::h() == hMinTest);
    REQUIRE(IceMinima::c() == cMinTest);
}
TEST_SUITE_END();
}
