/*!
 * @file eexample1-sasipmesh.cpp
 * @date 10 Jul 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "ParametricTools.hpp"
#include "SasipMesh.hpp"
#include "Interpolations.hpp"
#include "Mesh.hpp"
#include "ParametricTransport.hpp"
#include "dgLimiters.hpp"
#include "dgVisu.hpp"

#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "stopwatch.hpp"
#include "testtools.hpp"

namespace Nextsim {
extern Timer GlobalTimer;
}

bool WRITE_VTK = true; //!< set to true for vtk output
int WRITE_EVERY = 5;

double TOL = 1.e-10; //!< tolerance for checking test results

#define EDGEDOFS(DG) ((DG == 1) ? 1 : ((DG == 3) ? 2 : 3))

//! Packman-initial at 256000, 256000 with radius 128000
//! plus a smooth initial at 768000, 256000 with smaller radius
class PackmanPlus : public Nextsim::Interpolations::Function {

public:
    double operator()(double x, double y) const
    {
        //      return 1;
        double r = sqrt(pow((x - 240000.0) / 64000.0, 2.0) + pow((y - 120000.0) / 64000.0, 2.0));
        double r1 = sqrt(pow((x - 272000.0) / 96000.0, 2.0) + pow((y - 368000.0) / 96000.0, 2.0));
        if (r < 1.0)
            if ((16000.0 - fabs(x - 240000.0) < 0.0) || (y > 120000.0 - 16000.0))
                return 1.0 + exp(-5.0 * r1 * r1);
        return exp(-5.0 * r1 * r1);
    }
};

// Velocity
class InitialVX : public Nextsim::Interpolations::Function { // (0.5,0.2) m/s

    double _time;

public:
    void settime(double t) { _time = t; }

    double operator()(double x, double y) const
    {
        return (y / 256000 - 1.0) * 2.0 * M_PI;
    }
};
class InitialVY : public Nextsim::Interpolations::Function {

    double _time;

public:
    void settime(double t) { _time = t; }

    double operator()(double x, double y) const
    {
        return (1.0 - x / 256000.0) * 2.0 * M_PI;
    }
};

//////////////////////////////////////////////////

template <int DG>
class Test {

    const Nextsim::SasipMesh& smesh; //!< Stores a reference to the spacial mesh.

    size_t N; //!< size of mesh N x N

    size_t NT; //!< number of time steps
    double dt; //!< time step size

    Nextsim::Mesh mesh; //!< space mesh

    //! Velocity vectors and density
    Nextsim::CellVector<DG> phi, finalphi;

    //! Transport main class
    Nextsim::ParametricTransport<DG, EDGEDOFS(DG)> dgtransport;

    //! Velocity Field
    InitialVX VX;
    InitialVY VY;

    std::vector<double> values; //!< for storing numerical results (mass)
    std::vector<double> errors; //!< for storing error w.r.t. exact solution

    size_t writestep; //! write out n step in total (for debugging only)

public:
    Test(const Nextsim::SasipMesh& mesh)
        : smesh(mesh)
        , dgtransport(smesh)
        , writestep(40)
    {
        //! Set time stepping scheme. 2nd order for dg0 and dg1, 3rd order dG2
        if (DG < 3)
            dgtransport.settimesteppingscheme("rk2");
        else
            dgtransport.settimesteppingscheme("rk3");
    }

    Test() { std::cout << "call Test(N). N is number of mesh elements per row" << std::endl; }

    void init()
    {
        //! Init Time Mesh
        double cfl = 0.01;

        double hmin = smesh.hmin();
        double area = smesh.area();
        std::cout << "SasipMesh: hmin = " << hmin << "\tarea = " << area << std::endl;

        dt = cfl * hmin / 1.0; // max-velocity is 1
        double tmax = 256000.0;
        NT = (static_cast<int>((tmax / dt + 1) / 100 + 1) * 100); // No time steps dividable by 100
        std::cout << "dt = " << dt << "\t NT = " << NT << std::endl;
        //! Init Vectors
        phi.resize_by_mesh(smesh);
    }

    void run()
    {
        values.clear();

        Nextsim::GlobalTimer.reset();

        // init the test case, in particular resize vectors
        init();

        double area = 0.0;
        for (size_t i = 0; i < smesh.nelements; ++i)
            area += Nextsim::ParametricTools::massMatrix<DG>(smesh, i)(0, 0);
        std::cout << "Error in area: " << (area / 512000.0 / 512000.0 - 1.0) << std::endl;

        // initial density
        Nextsim::Interpolations::Function2DG(smesh, phi, PackmanPlus());

        if (WRITE_VTK) {
            Nextsim::VTK::write_dg<DG>("ResultsSasip/dg", 0, phi, smesh);
        }

        //! Save initial solution for error control
        finalphi = phi;

        //! time loop

        //	writestep = NT;
        for (size_t iter = 1; iter <= NT; ++iter) {
            VX.settime(dt * iter);
            VY.settime(dt * iter);
            Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVx(), VX);
	    Nextsim::Interpolations::Function2DG(smesh, dgtransport.GetVy(), VY);
            dgtransport.reinitnormalvelocity();

            dgtransport.step(dt, phi); // performs one time step with the 2nd Order Heun scheme
            if (WRITE_VTK)
                if (iter % (NT / writestep) == 0) {
                    Nextsim::VTK::write_dg<DG>("ResultsSasip/dg", iter / (NT / writestep), phi, smesh);
                    Nextsim::VTK::write_dg<DG>("ResultsSasip/vx", iter / (NT / writestep), dgtransport.GetVx(), smesh);
                    Nextsim::VTK::write_dg<DG>("ResultsSasip/vy", iter / (NT / writestep), dgtransport.GetVy(), smesh);
                }

            //	  if (iter==5) abort();
        }
        Nextsim::CellVector<DG> errorphi = phi;
        errorphi += -finalphi;
        if (WRITE_VTK) {
            Nextsim::VTK::write_dg<DG>("ResultsSasip/error", N, errorphi, smesh);
        }
        errors.push_back(errorphi.norm() * sqrt(mesh.hx * mesh.hy));
        values.push_back(phi.col(0).sum() * mesh.hx * mesh.hy);

        std::cout.precision(16);
    }

    void print_error(const std::string& message) const
    {
        std::cerr << "dG(" << DG << ") FAILED: " << message << std::endl;
        for (size_t i = 0; i < values.size(); ++i)
            std::cerr << values[i] << "\t" << errors[i] << std::endl;
    }

    void print_error(const std::array<double, 3>& v, const std::array<double, 3>& e,
        const std::string& message) const
    {
        std::cerr << "dG(" << DG << ") FAILED: " << message << std::endl;

        assert(values.size() >= 3);

        for (size_t i = 0; i < 3; ++i)
            std::cerr << v[i] << " = " << values[i + values.size() - 3] << "\t" << e[i] << " = "
                      << errors[i + values.size() - 3] << std::endl;
    }

    bool check_references(const std::array<double, 3>& v, const std::array<double, 3>& e) const
    {
        bool passed = true;

        for (size_t i = 0; i < 3; ++i) {
            if (fabs(v[i] - values[i]) > TOL) {
                print_error(v, e, "difference in mass");
                passed = false;
            }
            if (fabs(e[i] - errors[i]) > TOL) {
                print_error(v, e, "difference in error");
                passed = false;
            }
        }

        return passed;
    }

    bool check() const
    {
        std::array<double, 3> val_ref, err_ref;

        if (N == 80) {
            if (DG == 1) {
                val_ref = std::array<double, 3>(
                    { 0.1040973014386282, 0.1196042191148516, 0.1264035055795092 });
                err_ref = std::array<double, 3>(
                    { 0.2055757930605657, 0.1547718388297195, 0.1111151282733836 });
            } else if (DG == 3) {
                val_ref = std::array<double, 3>(
                    { 0.1297377725023032, 0.1290794989188777, 0.1290992303550586 });
                err_ref = std::array<double, 3>(
                    { 0.07138784362943089, 0.01941901084543833, 0.004572185826845114 });
            } else if (DG == 6) {
                val_ref = std::array<double, 3>(
                    { 0.1290847613765785, 0.129099700340886, 0.1290996998763931 });
                err_ref = std::array<double, 3>(
                    { 0.03113698402363593, 0.008527051921868931, 0.001234126542186417 });
            } else
                abort();
        } else {
            print_error("reference values only for N=20,40,80");
            return false;
        }

        return check_references(val_ref, err_ref);
    }
};

//////////////////////////////////////////////////

int main()
{
    Nextsim::SasipMesh smesh;
    smesh.readmesh("../SasipMesh/distortedrectangle.smesh");

    Test<3> test(smesh);

    test.run();
}
