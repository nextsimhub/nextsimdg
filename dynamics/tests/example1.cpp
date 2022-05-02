/*!
 * @file eexample1.cpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#include "Mesh.hpp"
#include "dgInitial.hpp"
#include "dgLimiters.hpp"
#include "dgTransport.hpp"
#include "dgVisu.hpp"

#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

//#include "stopwatch.hpp"
#include "testtools.hpp"

#include "../../src/include/Timer.hpp"

/*namespace Nextsim
{
  extern Timer GlobalTimer;
}*/

bool WRITE_VTK = false; //!< set to true for vtk output
int WRITE_EVERY = 5;

double TOL = 1.e-10; //!< tolerance for checking test results

#define EDGEDOFS(DG) ( (DG==1)?1:( (DG==3)?2:3) )

//! Initially a smooth bump centered at (0.4,0.4)
//! This will be transported in a circle with (-y, x) for one complete revolution
struct InitialPhi {

public:
    double smooth(double x) const // smooth transition from 0 to 1 on [0,1]
    {
        if (x <= 0)
            return 0;
        if (x >= 1.0)
            return 0;

        if (x < 0.5)
            return 0.5 * exp(-1. / x) / exp(-2.0);
        else
            return 1. - 0.5 * exp(-1. / (1. - x)) / exp(-2.0);
    }

    double operator()(double x, double y) const
    {
        double r = sqrt(pow(x - 0.4, 2.0) + pow(y - 0.4, 2.0));
        if (r < 0.1)
            return 1.0;
        if (r < 0.3)
            return 1.0 - smooth(5.0 * (r - 0.1));
        return 0.0;
    }
};

// Velocity
struct InitialVX {

    double _time;

public:
    void settime(double t) { _time = t; }

    double operator()(double x, double y) const
    {
        return (0.5 * M_PI * sin(0.5 * _time)) * ((y - 0.5));
    }
};
struct InitialVY {

    double _time;

public:
    void settime(double t) { _time = t; }

    double operator()(double x, double y) const
    {
        return (0.5 * M_PI * sin(0.5 * _time)) * (-(x - 0.5));
    }
};

//////////////////////////////////////////////////

template <int DG>
class Test {

    size_t N; //!< size of mesh N x N

    size_t NT; //!< number of time steps
    double dt; //!< time step size

    Nextsim::Mesh mesh; //!< space mesh

    //! Velocity vectors and density
    Nextsim::CellVector<DG> vx, vy, phi, finalphi;

    //! Transport main class
  Nextsim::DGTransport<DG,EDGEDOFS(DG)> dgtransport;

    //! Velocity Field
    InitialVX VX;
    InitialVY VY;

    std::vector<double> values; //!< for storing numerical results (mass)
    std::vector<double> errors; //!< for storing error w.r.t. exact solution

    size_t writestep; //! write out n step in total (for debugging only)

public:
    Test(size_t n)
        : N(n)
        , dgtransport(vx, vy)
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
        //! Init Mesh
        mesh.BasicInit(N, N, 1.0 / N, 1.0 / N);

        //! Init Time Mesh
        double cfl = 0.1;
        dt = cfl * std::min(mesh.hx, mesh.hy) / 1.0; // max-velocity is 1
        double tmax = 2.0 * M_PI;
        NT = (static_cast<int>((tmax / dt + 1) / 100 + 1) * 100); // No time steps dividable by 100

        //! Init Transport Scheme
        dgtransport.setmesh(mesh);

        //! Init Vectors
        vx.resize_by_mesh(mesh);
        vy.resize_by_mesh(mesh);
        phi.resize_by_mesh(mesh);
    }

    void run()
    {
        values.clear();

        size_t NITER = 3;

        //! Iteration with successive mesh refinement
        for (size_t iter = 0; iter < NITER; ++iter) {
            // Nextsim::GlobalTimer.reset();
            Nextsim::Timer::main.reset();
            // Nextsim::Timer::main.tick("run");
            Nextsim::Timer::main.tick("run");

            // Nextsim::Timer::main.tick("run -- init");
            Nextsim::Timer::main.tick("run -- init");

            init();

            // initial density
            Nextsim::L2ProjectInitial(mesh, phi, InitialPhi());

            if (WRITE_VTK) {
                Nextsim::Timer::main.tick("run -- init -- vtk");
                Nextsim::VTK::write_dg<DG>("Results/dg", 0, phi, mesh);
                Nextsim::Timer::main.tock("run -- init -- vtk");
            }

            //! Save initial solution for error control
            finalphi = phi;
            Nextsim::Timer::main.tock("run -- init");

            //! time loop
            Nextsim::Timer::main.tick("run -- loop");
            for (size_t iter = 1; iter <= NT; ++iter) {
                Nextsim::Timer::main.tick("run -- loop -- vel");
                VX.settime(dt * iter);
                VY.settime(dt * iter);
                Nextsim::Timer::main.tick("run -- loop -- vel -- l2");
                Nextsim::L2ProjectInitial(mesh, vx, VX);
                Nextsim::L2ProjectInitial(mesh, vy, VY);
                Nextsim::Timer::main.tock("run -- loop -- vel -- l2");
                dgtransport.reinitvelocity();
                Nextsim::Timer::main.tock("run -- loop -- vel");

                dgtransport.step(dt, phi); // performs one time step with the 2nd Order Heun scheme
                if (WRITE_VTK)
                    if (iter % (NT / writestep) == 0) {
                        Nextsim::Timer::main.tick("run -- loop -- vtk");
                        Nextsim::VTK::write_dg<DG>(
                            "Results/dg", iter / (NT / writestep), phi, mesh);
                        Nextsim::Timer::main.tock("run -- loop -- vtk");
                    }
            }
            Nextsim::Timer::main.tock("run -- loop");

            Nextsim::Timer::main.tick("run -- error");
            Nextsim::CellVector<DG> errorphi = phi;
            errorphi += -finalphi;
            if (WRITE_VTK) {
                Nextsim::Timer::main.tick("run -- error- vtk");
                Nextsim::VTK::write_dg<DG>("Results/error", N, errorphi, mesh);
                Nextsim::Timer::main.tock("run -- error- vtk");
            }
            errors.push_back(errorphi.norm() * sqrt(mesh.hx * mesh.hy));
            values.push_back(phi.col(0).sum() * mesh.hx * mesh.hy);
            Nextsim::Timer::main.tock("run -- error");

            Nextsim::Timer::main.tock("run");

            if (iter < NITER - 1)
                N *= 2; //!< double the problem size
        }

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
    size_t N = 20;

    Test<1> test0(N);
    test0.run();
    if (!test0.check())
        std::cout << "dG(0) TEST FAILED" << std::endl;
    else
        std::cout << "dG(0) TEST PASSED" << std::endl;

    Test<3> test1(N);
    test1.run();
    if (!test1.check())
        std::cout << "dG(1) TEST FAILED" << std::endl;
    else
        std::cout << "dG(1) TEST PASSED" << std::endl;

    Test<6> test2(N);
    test2.run();
    if (!test2.check())
        std::cout << "dG(2) TEST FAILED" << std::endl;
    else
        std::cout << "dG(2) TEST PASSED" << std::endl;

    std::cout << Nextsim::Timer::main << std::endl;
}
