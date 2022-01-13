#ifndef __BENCHMARK_DATA_HPP__
#define __BENCHMARK_DATA_HPP__

#include "dginitial.hpp"
#include "dynamics.hpp"

namespace ReferenceScale {
// Box-Test case from [Geosci. Model Dev., 8, 1747–1761, 2015]
constexpr double L = 1000000.0; //!< Size of domain
constexpr double vmax_ocean = 0.1; //!< Maximum velocity of ocean
constexpr double vmax_atm = 30.0 / exp(1.0); //!< Max. vel. of wind

constexpr double rho_ice = 900.0; //!< Sea ice density
constexpr double rho_atm = 1.3; //!< Air density
constexpr double rho_ocean = 1026.0; //!< Ocean density

constexpr double C_atm = 1.2e-3; //!< Air drag coefficient
constexpr double C_ocean = 5.5e-3; //!< Ocean drag coefficient

constexpr double Pstar = 27500; //!< Ice strength
constexpr double fc = 1.46e-4; //!< Coriolis
// parameters form nextsim options.cpp line 302
constexpr double compaction_param = -20; //!< Compation parameter
constexpr double undamaged_time_relaxation_sigma = 1e7; //!< seconds
constexpr double exponent_relaxation_sigma = 5;
constexpr double young = 5.9605e+08;
constexpr double nu0 = 1. / 3.; //!< \param Poisson's ratio
constexpr double compr_strength = 1e10; //! \param compr_strength (double) Maximum compressive strength [N/m2]
constexpr double tan_phi = 0.7; //! \param tan_phi (double) Internal friction coefficient (mu)
constexpr double C_lab = 2.0e6; //! \param C_lab (double) Cohesion at the lab scale (10 cm) [Pa]

}

inline constexpr double SQR(double x)
{
    return x * x;
}

//! Description of the problem data, wind & ocean fields
class OceanX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (2.0 * y / ReferenceScale::L - 1);
    }
};
class OceanY : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return ReferenceScale::vmax_ocean * (1.0 - 2.0 * x / ReferenceScale::L);
    }
};

class AtmX : virtual public Nextsim::InitialBase {
    double time;

public:
    void settime(double t)
    {
        time = t;
    }
    double operator()(double x, double y) const
    {
        double X = M_PI * x / ReferenceScale::L;
        double Y = M_PI * y / ReferenceScale::L;
        constexpr double T = 4.0 * 24.0 * 60.0 * 60.0; //!< 4 days
        return 5.0 + (sin(2 * M_PI * time / T) - 3.0) * sin(2 * X) * sin(Y);

        // return M_PI * M_PI / 2.0 / SQR(ReferenceScale::L)
        //     * (3.0 * sin(X) * sin(Y) - 4.0 * cos(2 * X) * cos(2 * Y));

        //! Center of cyclone (in m)
        double cM = 256000. + 51200. * time / (24.0 * 60.0 * 60.0);

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cM) + SQR(y - cM))) * 1.e-3;

        double alpha = 72.0 / 180.0 * M_PI;
        return -scale * ReferenceScale::vmax_atm * (cos(alpha) * (x - cM) + sin(alpha) * (y - cM));
    }
};
class AtmY : virtual public Nextsim::InitialBase {
    double time;

public:
    void settime(double t)
    {
        time = t;
    }
    double operator()(double x, double y) const
    {
        double X = M_PI * x / ReferenceScale::L;
        double Y = M_PI * y / ReferenceScale::L;
        constexpr double T = 4.0 * 24.0 * 60.0 * 60.0; //!< 4 days
        return 5.0 + (sin(2 * M_PI * time / T) - 3) * sin(2 * Y) * sin(X);

        // return M_PI * M_PI / 2.0 / SQR(ReferenceScale::L)
        //     * (-cos(X) * cos(Y) + 12.0 * sin(2 * X) * sin(2 * Y));

        //! Center of cyclone (in m)
        double cM = 256000. + 51200. * time / (24.0 * 60.0 * 60.0);

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01e-3 * sqrt(SQR(x - cM) + SQR(y - cM))) * 1.e-3;

        double alpha = 72.0 / 180.0 * M_PI;
        return -scale * ReferenceScale::vmax_atm * (-sin(alpha) * (x - cM) + cos(alpha) * (y - cM));
    }
};
class ExX : virtual public Nextsim::InitialBase {

public:
    double operator()(double x, double y) const
    {
        double X = M_PI * x / ReferenceScale::L;
        double Y = M_PI * y / ReferenceScale::L;
        return sin(X) * sin(Y);
    }
};
class ExY : virtual public Nextsim::InitialBase {

public:
    double operator()(double x, double y) const
    {
        double X = 2.0 * M_PI * x / ReferenceScale::L;
        double Y = 2.0 * M_PI * y / ReferenceScale::L;
        return sin(X) * sin(Y);
    }
};

class InitialH : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        // x and y are given in meters
        return 0.3 + 0.005 * (sin(6.e-5 * x) + sin(3.e-5 * y));
    }
};
class InitialA : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return x / ReferenceScale::L;

        return 1.0;
    }
};

#endif
