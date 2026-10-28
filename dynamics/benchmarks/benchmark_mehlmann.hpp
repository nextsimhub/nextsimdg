/*!
 * @file benchmark_mehlmann.hpp
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 * The sea ice benchmark of Mehlmann et al. (2021), "Simulating linear kinematic features in
 * viscous-plastic sea ice models on quadrilateral and triangular grids with different variable
 * staggering", JAMES. Driver for standalone executables that run a dynamics kernel
 * (CG mEVP, BBM, LDG mEVP) without the rest of the model.
 *
 * Domain [0, 512 km]^2 with one layer of land cells around it (as run/make_init_benchmark.py),
 * ice concentration 1, ice height 0.3 + 0.005 (sin(6e-5 x) + sin(3e-5 y)), circular ocean
 * current and a cyclone moving to the north east, as in the modules BenchmarkOcean and
 * BenchmarkAtmosphere of the full model. Forcing and velocities are cell averages (dG0), as in
 * the model.
 */

#ifndef BENCHMARK_MEHLMANN_HPP
#define BENCHMARK_MEHLMANN_HPP

#include "DGModelArray.hpp"
#include "Interpolations.hpp"
#include "ParametricMesh.hpp"
#include "dgVector.hpp"
#include "dgVisu.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"
#include "include/gridNames.hpp"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace Nextsim {
namespace Benchmark {

    constexpr FloatType L = 512000.0; //!< size of the domain [m]
    constexpr FloatType oneDay = 24.0 * 60.0 * 60.0;
    constexpr FloatType vMaxOcean = 0.01; //!< maximum velocity of the ocean [m/s]
    constexpr FloatType vMaxAtm = 30.0; //!< maximum velocity of the wind [m/s]
    constexpr FloatType cycloneDuration = 5.0; //!< days for the cyclone to cross the domain

    //! Settings of a benchmark run, set from the command line
    struct Settings {
        size_t nx = 64; //!< cells in each direction
        FloatType dt = 120.0; //!< time step [s]
        FloatType days = 2.0; //!< simulated time
        FloatType outputHours = 6.0; //!< interval of the VTK output, 0: no output
        int nSteps = 120; //!< mEVP / BBM sub-iterations per time step
        std::string directory; //!< output directory (default: name of the benchmark)
    };

    inline void usage(const std::string& name)
    {
        std::cout << "Usage: " << name
                  << " [--nx N] [--dt seconds] [--days d] [--output hours] [--nsteps n] [--dir "
                     "directory]"
                  << std::endl;
    }

    inline Settings parseArguments(int argc, char* argv[], const std::string& name)
    {
        Settings s;
        s.directory = name;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                usage(argv[0]);
                std::exit(0);
            }
            if (i + 1 >= argc) {
                usage(argv[0]);
                throw std::invalid_argument("missing value for " + arg);
            }
            const std::string value = argv[++i];
            if (arg == "--nx")
                s.nx = std::stoul(value);
            else if (arg == "--dt")
                s.dt = std::stod(value);
            else if (arg == "--days")
                s.days = std::stod(value);
            else if (arg == "--output")
                s.outputHours = std::stod(value);
            else if (arg == "--nsteps")
                s.nSteps = std::stoi(value);
            else if (arg == "--dir")
                s.directory = value;
            else {
                usage(argv[0]);
                throw std::invalid_argument("unknown argument " + arg);
            }
        }
        return s;
    }

    //! Initial ice height
    class InitialH : public Interpolations::Function {
    public:
        FloatType operator()(FloatType x, FloatType y) const
        {
            return 0.3 + 0.005 * (std::sin(6.e-5 * x) + std::sin(3.e-5 * y));
        }
    };
    //! Constant function
    class Constant : public Interpolations::Function {
        FloatType value;

    public:
        Constant(FloatType v)
            : value(v)
        {
        }
        FloatType operator()(FloatType x, FloatType y) const { return value; }
    };

    /*!
     * The data of the benchmark as ModelArrays, as the model passes them to the dynamics
     * kernels. Velocities and forcing are cell averages, ice fields are DG fields.
     */
    class Problem {
    public:
        const Settings settings;
        ParametricMesh smesh; //!< the mesh, for the initial values and the output
        ModelArray coords, mask;
        ModelArray hice, cice, hsnow, damage; //!< DG fields
        ModelArray uWind, vWind, uOcean, vOcean, ssh, u0, v0;

        Problem(const Settings& s)
            : settings(s)
            , smesh(CARTESIAN)
        {
            const size_t nx = s.nx;
            ModelArray::setDimension(ModelArray::Dimension::X, nx);
            ModelArray::setDimension(ModelArray::Dimension::Y, nx);
            ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
            ModelArray::setDimension(ModelArray::Dimension::YVERTEX, nx + 1);
            ModelArray::setDimension(ModelArray::Dimension::XCG, CGDEGREE * nx + 1);
            ModelArray::setDimension(ModelArray::Dimension::YCG, CGDEGREE * nx + 1);

            coords = ModelArray(ModelArray::Type::VERTEX);
            coords.reinitialize();
            for (size_t iy = 0, idx = 0; iy <= nx; ++iy)
                for (size_t ix = 0; ix <= nx; ++ix, ++idx) {
                    coords.components(idx)[0] = L * ix / nx;
                    coords.components(idx)[1] = L * iy / nx;
                }

            // ice everywhere, one layer of land around the domain
            mask = ModelArray(ModelArray::Type::H);
            mask.reinitialize();
            for (size_t iy = 0, idx = 0; iy < nx; ++iy)
                for (size_t ix = 0; ix < nx; ++ix, ++idx)
                    mask[idx] = (ix == 0 || iy == 0 || ix + 1 == nx || iy + 1 == nx) ? 0.0 : 1.0;

            smesh.coordinatesFromModelArray(coords);
            smesh.landmaskFromModelArray(mask);

            // initial ice fields, L2 projection onto the DG space
            hice = projectDG(InitialH(), true);
            cice = projectDG(Constant(1.0), true);
            hsnow = projectDG(Constant(0.0), true);
            damage = projectDG(Constant(1.0), false);

            for (ModelArray* f : { &uWind, &uOcean, &u0 })
                *f = ModelArray(ModelArray::Type::U);
            for (ModelArray* f : { &vWind, &vOcean, &v0 })
                *f = ModelArray(ModelArray::Type::V);
            ssh = ModelArray(ModelArray::Type::H);
            for (ModelArray* f : { &uWind, &vWind, &uOcean, &vOcean, &ssh, &u0, &v0 }) {
                f->reinitialize();
                *f = 0.0;
            }

            // circular ocean current, constant in time
            for (size_t iy = 0, idx = 0; iy < nx; ++iy)
                for (size_t ix = 0; ix < nx; ++ix, ++idx) {
                    const FloatType fx = (ix + 0.5) / nx;
                    const FloatType fy = (iy + 0.5) / nx;
                    uOcean[idx] = mask[idx] * vMaxOcean * (2.0 * fy - 1.0);
                    vOcean[idx] = mask[idx] * vMaxOcean * (1.0 - 2.0 * fx);
                }
        }

        //! The cyclone at time t [s], in the cell centres
        void setWind(FloatType t)
        {
            const size_t nx = settings.nx;
            const FloatType A = 1.e-5;
            const FloatType k = 1.e-5; // 1/m
            const FloatType alpha = 72. / 180. * M_PI;
            // centre of the cyclone
            const FloatType c = 0.5 * L * (1.0 + t / oneDay / cycloneDuration);
            for (size_t iy = 0, idx = 0; iy < nx; ++iy)
                for (size_t ix = 0; ix < nx; ++ix, ++idx) {
                    const FloatType xp = L * (ix + 0.5) / nx - c;
                    const FloatType yp = L * (iy + 0.5) / nx - c;
                    const FloatType scale = A * std::exp(-k * std::hypot(xp, yp));
                    uWind[idx] = -scale * vMaxAtm * (std::cos(alpha) * xp + std::sin(alpha) * yp);
                    vWind[idx] = -scale * vMaxAtm * (-std::sin(alpha) * xp + std::cos(alpha) * yp);
                }
        }

        //! DG field as DGVector
        DGVector<DGCOMP> dg(const ModelArray& field) const
        {
            DGVector<DGCOMP> v(smesh);
            DGModelArray::ma2dg(field, v);
            return v;
        }

        //! Cell averages as DGVector<1>
        DGVector<1> dg0(const ModelArray& field) const
        {
            DGVector<1> v(smesh);
            DGModelArray::ma2dg(field, v);
            return v;
        }

        //! integral of the ice height (DG0 component times area)
        FloatType iceVolume() const
        {
            const DGVector<DGCOMP> h = dg(hice);
            FloatType vol = 0;
            for (size_t i = 0; i < smesh.nelements; ++i)
                vol += h(i, 0) * smesh.area(i);
            return vol;
        }

    private:
        //! L2 projection of a function onto the DG space, set to 0 on land
        ModelArray projectDG(const Interpolations::Function& f, bool masked)
        {
            DGVector<DGCOMP> v(smesh);
            Interpolations::Function2DG(smesh, v, f);
            if (masked)
                for (size_t i = 0; i < smesh.nelements; ++i)
                    if (!smesh.landmask[i])
                        v.row(i).setZero();
            ModelArray ma(ModelArray::Type::DG);
            ma.reinitialize();
            return DGModelArray::dg2ma(v, ma);
        }
    };

    /*!
     * Runs the benchmark with a dynamics kernel. withDamage: the kernel uses the damage field
     * (BBM).
     */
    template <class Kernel>
    void run(Kernel& kernel, Problem& problem, const std::string& name, bool withDamage = false)
    {
        const Settings& s = problem.settings;
        const size_t nt = std::lround(s.days * oneDay / s.dt);
        const size_t ntOutput = (s.outputHours > 0)
            ? std::max<long>(1, std::lround(s.outputHours * 3600. / s.dt))
            : 0;
        const size_t ntLog = std::max<long>(1, std::lround(3600. / s.dt));

        std::cout << name << ": Mehlmann benchmark, " << s.nx << " x " << s.nx
                  << " cells, dt = " << s.dt << " s, " << nt << " time steps, " << s.nSteps
                  << " sub-iterations, advection DG" << DGCOMP << std::endl;

        kernel.initialise(problem.coords, false, problem.mask);
        kernel.setData(uName, problem.u0);
        kernel.setData(vName, problem.v0);

        if (ntOutput > 0)
            std::filesystem::create_directories(s.directory);
        auto output = [&](size_t step) {
            if (ntOutput == 0)
                return;
            const int n = step / ntOutput;
            const std::string d = s.directory + "/";
            VTK::write_dg(d + "H", n, problem.dg(problem.hice), problem.smesh);
            VTK::write_dg(d + "A", n, problem.dg(problem.cice), problem.smesh);
            if (withDamage)
                VTK::write_dg(d + "D", n, problem.dg(problem.damage), problem.smesh);
            VTK::write_dg(d + "u", n, problem.dg0(kernel.getDG0Data(uName)), problem.smesh);
            VTK::write_dg(d + "v", n, problem.dg0(kernel.getDG0Data(vName)), problem.smesh);
            VTK::write_dg(d + "shear", n, problem.dg0(kernel.getDG0Data(shearName)), problem.smesh);
        };

        const FloatType volume0 = problem.iceVolume();
        const TimePoint start("2023-01-01T00:00:00Z");
        const Duration step(s.dt);
        const auto wallStart = std::chrono::steady_clock::now();

        output(0);
        for (size_t it = 0; it < nt; ++it) {
            const FloatType t = it * s.dt;
            TimePoint now = start;
            now += Duration(t);
            const TimestepTime tst = { now, step };

            // the order of the calls as in the model (PrognosticData::update)
            problem.setWind(t);
            kernel.setDGArray(hiceName, static_cast<ModelArray::DataType&>(problem.hice));
            kernel.setDGArray(ciceName, static_cast<ModelArray::DataType&>(problem.cice));
            kernel.setDGArray(hsnowName, static_cast<ModelArray::DataType&>(problem.hsnow));
            if (withDamage)
                kernel.setDGArray(damageName, static_cast<ModelArray::DataType&>(problem.damage));
            kernel.setData(uWindName, problem.uWind);
            kernel.setData(vWindName, problem.vWind);
            kernel.setData(uOceanName, problem.uOcean);
            kernel.setData(vOceanName, problem.vOcean);
            kernel.setData(sshName, problem.ssh);
            kernel.prepareAdvection();
            kernel.update(tst);

            if ((it + 1) % ntLog == 0 || it + 1 == nt) {
                const DGVector<1> u = problem.dg0(kernel.getDG0Data(uName));
                const DGVector<1> v = problem.dg0(kernel.getDG0Data(vName));
                const FloatType umax = (u.array().square() + v.array().square()).sqrt().maxCoeff();
                const FloatType seconds
                    = std::chrono::duration<FloatType>(std::chrono::steady_clock::now() - wallStart)
                          .count();
                std::cout << std::fixed << std::setprecision(2) << "t = " << std::setw(6)
                          << (it + 1) * s.dt / 3600. << " h  |u|max = " << std::setprecision(4)
                          << umax << " m/s  volume change = " << std::scientific
                          << std::setprecision(3) << problem.iceVolume() / volume0 - 1.0
                          << std::fixed << std::setprecision(1) << "  wall " << seconds << " s"
                          << std::endl;
            }
            if (ntOutput > 0 && (it + 1) % ntOutput == 0)
                output(it + 1);
        }
    }

} // namespace Benchmark
} // namespace Nextsim

#endif /* BENCHMARK_MEHLMANN_HPP */
