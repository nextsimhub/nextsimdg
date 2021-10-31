#ifndef __BENCHMARK_DATA_HPP__
#define __BENCHMARK_DATA_HPP__

#include "dginitial.hpp"
#include "dynamics.hpp"

inline double SQR(double x) 
{
  return x*x;
}


//! Reference scaling 
class ReferenceScale
{
public:
  static constexpr double L = 512.e3; // 1  |-> 512 km
  static constexpr double T = 1.e3;   // 1  |-> 1000 sec approx 16 min.
};


//! Description of the problem data, wind & ocean fields
class OceanX : virtual public Nextsim::InitialBase {
public:
  double operator()(double x, double y) const
  {
    double Y = y*ReferenceScale::L; //!< Coordinate in m

    //! maximum velocity in reference system
    double vmax = 0.01 * ReferenceScale::T/ReferenceScale::L;
    return vmax * (2.0*Y/ReferenceScale::L-1);
  }
};
class OceanY : virtual public Nextsim::InitialBase {
public:
  double operator()(double x, double y) const
  {
    double X = x*ReferenceScale::L; //!< Coordinate in m

    //! maximum velocity in reference system
    double vmax = 0.01 * ReferenceScale::T/ReferenceScale::L;
    return vmax * (1.0 - 2.0*X/ReferenceScale::L);
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
    //! Center of cyclone (in km) 
    double cKM = 256.0 + 51.2 * time * ReferenceScale::T / (24.0*60.0*60.0);

    //! coordinate (in km)
    double xKM = x * ReferenceScale::L * 1.e-3;
    double yKM = y * ReferenceScale::L * 1.e-3;

    //! maximum velocity (in reference system)
    double vmax = 30.0/exp(1.0) * ReferenceScale::T/ReferenceScale::L;

    //! scaling factor to reduce wind away from center
    double scale = exp(1.0)/100.0 * exp(-0.01 * sqrt( SQR(xKM-cKM) + SQR(yKM-cKM) ));

    double alpha = 72.0/180.0*M_PI;
    return -scale * vmax * (cos(alpha)*(xKM-cKM) + sin(alpha)*(yKM-cKM));
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
    //! Center of cyclone (in km) 
    double cKM = 256.0 + 51.2 * time * ReferenceScale::T / (24.0*60.0*60.0);

    //! coordinate (in km)
    double xKM = x * ReferenceScale::L * 1.e-3;
    double yKM = y * ReferenceScale::L * 1.e-3;

    //! maximum velocity (in reference system)
    double vmax = 30.0/exp(1.0) * ReferenceScale::T/ReferenceScale::L;

    //! scaling factor to reduce wind away from center
    double scale = exp(1.0)/100.0 * exp(-0.01 * sqrt( SQR(xKM-cKM) + SQR(yKM-cKM) ));

    double alpha = 72.0/180.0*M_PI;
    return -scale * vmax * (-sin(alpha)*(xKM-cKM) + cos(alpha)*(yKM-cKM));
  }
};


class InitialH : virtual public Nextsim::InitialBase
{
public:
  double operator()(double x, double y) const
  {
    //! coordinate (in km)
    double xKM = x * ReferenceScale::L * 1.e-3;
    double yKM = y * ReferenceScale::L * 1.e-3;
    return 0.3 + 0.005 * (sin(1.e-3 * 60.0 * xKM)+sin(1.e-3 * 30.0 * yKM));
  }
};
class InitialA : virtual public Nextsim::InitialBase
{
public:
  double operator()(double x, double y) const
  {
    return 1.0;
  }
};

#endif

