/*!
 * @file BenchmarkAtmosphere.cpp
 *
 * @date 24 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BenchmarkAtmosphere.hpp"

#include "include/BenchmarkCoordinates.hpp"
#include "include/Time.hpp"

namespace Nextsim {

void BenchmarkAtmosphere::setData(const ModelState::DataMap& ms)
{
    IAtmosphereBoundary::setData(ms);
    BenchmarkCoordinates::setData();
    // Constant, zero fluxes in the atmosphere
    qia = 0.;
    dqia_dt = 0.;
    qow = 0.;
    subl = 0.;
    snow = 0.;
    rain = 0.;
    evap = 0.;
}

void BenchmarkAtmosphere::update(const TimestepTime& tst)
{
    IAtmosphereBoundary::update(tst);

    // set the initial time on the first update
    if (!t0Set) {
        t0 = tst.start;
        t0Set = true;
    }
    // length of 1 day in seconds
    constexpr double oneday = 24.0 * 60.0 * 60.0;
    // maximum wind velocity of the cyclone
    constexpr double vMax = 30.0;

    // number of days elapsed since t0
    Duration elapsedTime = tst.start - t0;
    double timeFraction = elapsedTime.seconds() / oneday;
    double cycloneDuration = 5.; // days

    // cyclone parameters
    const double A = 1e-5;
    const double k = 1e-5; // m⁻¹
    const double alpha = 72. / 180. * M_PI;
    const double cosalpha = cos(alpha);
    const double sinalpha = sin(alpha);

    // centre of the cyclone (metres)
    double x0 = BenchmarkCoordinates::nx() * BenchmarkCoordinates::dx * 0.5
        * (1 + timeFraction / cycloneDuration);
    double y0 = BenchmarkCoordinates::ny() * BenchmarkCoordinates::dy * 0.5
        * (1 + timeFraction / cycloneDuration);

    // distance from the centre of the cyclone
    const ModelArray& xPrime = BenchmarkCoordinates::x() - x0;
    const ModelArray& yPrime = BenchmarkCoordinates::y() - y0;

    // Perform the rest of the calculation per element
    for (size_t j = 0; j < BenchmarkCoordinates::ny(); ++j) {
        for (size_t i = 0; i < BenchmarkCoordinates::nx(); ++i) {
            // Expression taken from the original implementation:
            // double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cMx) + SQR(y - cMy)))
            // * 1.e-3;
            double scale = A * exp(-k * hypot(xPrime(i, j), yPrime(i, j)));
            uwind(i, j) = -scale * vMax * (cosalpha * xPrime(i, j) + sinalpha * yPrime(i, j));
            vwind(i, j) = -scale * vMax * (-sinalpha * xPrime(i, j) + cosalpha * yPrime(i, j));
        }
    }
}
} /* namespace Nextsim */
