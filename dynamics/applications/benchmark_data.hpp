#ifndef __BENCHMARK_DATA_HPP__
#define __BENCHMARK_DATA_HPP__

#include "dginitial.hpp"
#include "dynamics.hpp"
#include "referencescale.hpp"

inline double SQR(double x)
{
    return x * x;
}

//! Description of the problem data, wind & ocean fields
class OceanX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 1.; //y*y;
        double Y = y * Nextsim::ReferenceScale::L; //!< Coordinate in m

        //! maximum velocity in reference system
        double vmax = 0.01 * Nextsim::ReferenceScale::T / Nextsim::ReferenceScale::L;
        return vmax * (2.0 * Y / Nextsim::ReferenceScale::L - 1);
    }
};
class OceanY : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 1; //x*x;
        double X = x * Nextsim::ReferenceScale::L; //!< Coordinate in m

        //! maximum velocity in reference system
        double vmax = 0.01 * Nextsim::ReferenceScale::T / Nextsim::ReferenceScale::L;
        return vmax * (1.0 - 2.0 * X / Nextsim::ReferenceScale::L);
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
        double sx = sin(M_PI * x);
        double sy = sin(M_PI * y);
        double s2x = sin(2. * M_PI * x);
        double s2y = sin(2. * M_PI * y);
        double cx = cos(M_PI * x);
        double cy = cos(M_PI * y);
        double c2x = cos(2. * M_PI * x);
        double c2y = cos(2. * M_PI * y);
        return 3.0 * M_PI * M_PI / 2. * sx * sy - 2.0 * M_PI * M_PI * c2x * c2y;

        return 1.0;

        return 1.e-5;
        //! Center of cyclone (in km)
        double cKM = 256.0 + 51.2 * time * Nextsim::ReferenceScale::T / (24.0 * 60.0 * 60.0);

        //! coordinate (in km)
        double xKM = x * Nextsim::ReferenceScale::L * 1.e-3;
        double yKM = y * Nextsim::ReferenceScale::L * 1.e-3;

        //! maximum velocity (in reference system)
        double vmax = 30.0 / exp(1.0) * Nextsim::ReferenceScale::T / Nextsim::ReferenceScale::L;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01 * sqrt(SQR(xKM - cKM) + SQR(yKM - cKM)));

        double alpha = 72.0 / 180.0 * M_PI;
        return -scale * vmax * (cos(alpha) * (xKM - cKM) + sin(alpha) * (yKM - cKM));
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

        double sx = sin(M_PI * x);
        double sy = sin(M_PI * y);
        double s2x = sin(2. * M_PI * x);
        double s2y = sin(2. * M_PI * y);
        double cx = cos(M_PI * x);
        double cy = cos(M_PI * y);
        double c2x = cos(2. * M_PI * x);
        double c2y = cos(2. * M_PI * y);
        return -M_PI * M_PI / 2. * cx * cy + 6.0 * M_PI * M_PI * s2x * s2y;

        return 1.0;

        return 1.e-5;
        //! Center of cyclone (in km)
        double cKM = 256.0 + 51.2 * time * Nextsim::ReferenceScale::T / (24.0 * 60.0 * 60.0);

        //! coordinate (in km)
        double xKM = x * Nextsim::ReferenceScale::L * 1.e-3;
        double yKM = y * Nextsim::ReferenceScale::L * 1.e-3;

        //! maximum velocity (in reference system)
        double vmax = 30.0 / exp(1.0) * Nextsim::ReferenceScale::T / Nextsim::ReferenceScale::L;

        //! scaling factor to reduce wind away from center
        double scale = exp(1.0) / 100.0 * exp(-0.01 * sqrt(SQR(xKM - cKM) + SQR(yKM - cKM)));

        double alpha = 72.0 / 180.0 * M_PI;
        return -scale * vmax * (-sin(alpha) * (xKM - cKM) + cos(alpha) * (yKM - cKM));
    }
};

class InitialVX : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return sin(M_PI * x) * sin(M_PI * y) + sin(M_PI * 10.0 * x) * sin(M_PI * 5 * y) * 0.3;
    }
};
class InitialVY : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y) + sin(M_PI * 3. * x) * sin(M_PI * 7 * y) * 0.4;
    }
};

class InitialH : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        //! coordinate (in km)
        double xKM = x * Nextsim::ReferenceScale::L * 1.e-3;
        double yKM = y * Nextsim::ReferenceScale::L * 1.e-3;
        return 0.3 + 0.005 * (sin(1.e-3 * 60.0 * xKM) + sin(1.e-3 * 30.0 * yKM));
    }
};
class InitialA : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 1.0;
    }
};

class InitialS11 : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    { //! coordinate (in km)
        double xKM = x * Nextsim::ReferenceScale::L * 1.e-3;
        double yKM = y * Nextsim::ReferenceScale::L * 1.e-3;
        return 0.0; // - 1e-1*sqrt(xKM*xKM + yKM*yKM);
    }
};
class InitialS12 : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    {
        return 0.0;
    }
};
class InitialS22 : virtual public Nextsim::InitialBase {
public:
    double operator()(double x, double y) const
    { //! coordinate (in km)
        double xKM = x * Nextsim::ReferenceScale::L * 1.e-3;
        double yKM = y * Nextsim::ReferenceScale::L * 1.e-3;
        return 0.0; // + 1e-1*sqrt(xKM*xKM + yKM*yKM) ;
    }
};

#endif
