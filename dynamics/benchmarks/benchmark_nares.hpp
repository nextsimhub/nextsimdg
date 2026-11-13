/*!
 * @file benchmark_nares.hpp
 *
 * Ice channel (Nares Strait) test case of Mehlmann et al. (2021). Driver for standalone
 * executables that run a dynamics kernel (CG mEVP, LDG mEVP) without the rest of the model.
 *
 * Same set-up as the full model with run/config_nares.cfg and run/make_init_nares.py:
 * channel [0, 240 km] x [0, 120 km], land in the first and last row of cells and a narrowing
 * in the middle, the left and right ends of the channel are open (no land cells there).
 * Ice concentration 1, ice height 0.3 + 0.005 (sin(6e-5 x) + sin(3e-5 y)) as in the box
 * benchmark, ocean at rest (ConstantOceanBoundary). Wind: unlike run/config_nares.cfg (cyclone
 * of BenchmarkAtmosphere), a uniform wind along the channel from left to right (+x), increased
 * linearly from 0 to windSpeed during windRamp. Forcing and velocities are cell averages (dG0),
 * as in the model.
 *
 * The command line options and the initial ice height come from benchmark_mehlmann.hpp. --nx N is the number of cells across the channel (y direction), the
 * channel has 2 N cells along it (x direction). Default N = 64: 128 x 64 cells, as
 * run/config_nares.cfg.
 */

#ifndef BENCHMARK_NARES_HPP
#define BENCHMARK_NARES_HPP

#include "benchmark_mehlmann.hpp"

#include <algorithm>

namespace Nextsim {
namespace Nares {

    constexpr FloatType Lx = 240000.0; //!< length of the channel [m]
    constexpr FloatType Ly = 120000.0; //!< width of the channel [m]
    constexpr FloatType windSpeed = 10.0; //!< wind along the channel, left to right [m/s]
    constexpr FloatType windRamp = 6.0 * 3600.0; //!< time for the wind to reach windSpeed [s]

    /*!
     * The data of the channel as ModelArrays, as the model passes them to the dynamics
     * kernels. Velocities and forcing are cell averages, ice fields are DG fields.
     */
    class Problem {
    public:
        const Benchmark::Settings settings;
        const size_t nx, ny; //!< cells along (x) and across (y) the channel
        ParametricMesh smesh; //!< the mesh, for the initial values and the output
        ModelArray coords, mask;
        ModelArray hice, cice, hsnow; //!< DG fields
        ModelArray uWind, vWind, uOcean, vOcean, ssh, u0, v0;

        Problem(const Benchmark::Settings& s)
            : settings(s)
            , nx(2 * s.nx)
            , ny(s.nx)
            , smesh(CARTESIAN)
        {
            ModelArray::setDimension(ModelArray::Dimension::X, nx);
            ModelArray::setDimension(ModelArray::Dimension::Y, ny);
            ModelArray::setDimension(ModelArray::Dimension::XVERTEX, nx + 1);
            ModelArray::setDimension(ModelArray::Dimension::YVERTEX, ny + 1);
            ModelArray::setDimension(ModelArray::Dimension::XCG, CGDEGREE * nx + 1);
            ModelArray::setDimension(ModelArray::Dimension::YCG, CGDEGREE * ny + 1);

            coords = ModelArray(ModelArray::Type::VERTEX);
            coords.reinitialize();
            for (size_t iy = 0, idx = 0; iy <= ny; ++iy)
                for (size_t ix = 0; ix <= nx; ++ix, ++idx) {
                    coords.components(idx)[0] = Lx * ix / nx;
                    coords.components(idx)[1] = Ly * iy / ny;
                }

            mask = ModelArray(ModelArray::Type::H);
            mask.reinitialize();
            for (size_t iy = 0, idx = 0; iy < ny; ++iy)
                for (size_t ix = 0; ix < nx; ++ix, ++idx)
                    mask[idx] = isIce(ix, iy) ? 1.0 : 0.0;

            smesh.coordinatesFromModelArray(coords);
            smesh.landmaskFromModelArray(mask);

            // initial ice fields, L2 projection onto the DG space
            hice = projectDG(Benchmark::InitialH());
            cice = projectDG(Benchmark::Constant(1.0));
            hsnow = projectDG(Benchmark::Constant(0.0));

            for (ModelArray* f : { &uWind, &uOcean, &u0 })
                *f = ModelArray(ModelArray::Type::U);
            for (ModelArray* f : { &vWind, &vOcean, &v0 })
                *f = ModelArray(ModelArray::Type::V);
            ssh = ModelArray(ModelArray::Type::H);
            // ocean at rest, ice at rest at t = 0
            for (ModelArray* f : { &uWind, &vWind, &uOcean, &vOcean, &ssh, &u0, &v0 }) {
                f->reinitialize();
                *f = 0.0;
            }
        }

        /*!
         * Uniform wind along the channel, from left to right (+x), at time t [s]:
         * u = windSpeed * min(t / windRamp, 1), v = 0.
         */
        void setWind(FloatType t)
        {
            const FloatType ramp = (windRamp > 0) ? std::min(t / windRamp, 1.0) : 1.0;
            uWind = windSpeed * ramp;
            vWind = 0.0;
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
        /*!
         * Mask of run/make_init_nares.py (lengths in km): land in the first and last row of
         * cells, left and right ends open, and the narrowing
         *
         *  (0,120)---- (80,120)                      / (240,120)
         *                       \                 /
         *                       (130,80) -(170,80)
         *
         *                       (130,40) -(170,40)
         *                       /                 \
         *  (0,0)------ (80,0)                       \  (240,0)
         *
         * land if dy < 20, x > 80 + 5/4 dy and x < 170 + 7/4 (40 - dy), with dy the distance
         * to the nearest wall y = 0 or y = 120, in the cell centres.
         */
        bool isIce(size_t ix, size_t iy) const
        {
            if (iy == 0 || iy + 1 == ny)
                return false;
            const FloatType x = 1.e-3 * Lx * (ix + 0.5) / nx; // km
            const FloatType y = 1.e-3 * Ly * (iy + 0.5) / ny; // km
            const FloatType dy = std::min(y, 1.e-3 * Ly - y);
            const bool land
                = dy < 20.0 && x > 80.0 + 5.0 / 4.0 * dy && x < 170.0 + 7.0 / 4.0 * (40.0 - dy);
            return !land;
        }

        //! L2 projection of a function onto the DG space, set to 0 on land
        ModelArray projectDG(const Interpolations::Function& f)
        {
            DGVector<DGCOMP> v(smesh);
            Interpolations::Function2DG(smesh, v, f);
            for (size_t i = 0; i < smesh.nelements; ++i)
                if (!smesh.landmask[i])
                    v.row(i).setZero();
            ModelArray ma(ModelArray::Type::DG);
            ma.reinitialize();
            return DGModelArray::dg2ma(v, ma);
        }
    };

    /*!
     * Runs the channel with a dynamics kernel, time loop as Benchmark::run (without damage).
     * VTK output: H, A (DG fields), u, v, speed V = |u| and shear (cell averages).
     */
    template <class Kernel> void run(Kernel& kernel, Problem& problem, const std::string& name)
    {
        const Benchmark::Settings& s = problem.settings;
        const size_t nt = std::lround(s.days * Benchmark::oneDay / s.dt);
        const size_t ntOutput = (s.outputHours > 0)
            ? std::max<long>(1, std::lround(s.outputHours * 3600. / s.dt))
            : 0;
        const size_t ntLog = std::max<long>(1, std::lround(3600. / s.dt));

        std::cout << name << ": Nares channel, " << problem.nx << " x " << problem.ny
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
            const DGVector<1> u = problem.dg0(kernel.getDG0Data(uName));
            const DGVector<1> v = problem.dg0(kernel.getDG0Data(vName));
            // speed |u| = sqrt(u^2 + v^2), cell averages
            DGVector<1> speed(problem.smesh);
            speed.array() = (u.array().square() + v.array().square()).sqrt();
            VTK::write_dg(d + "u", n, u, problem.smesh);
            VTK::write_dg(d + "v", n, v, problem.smesh);
            VTK::write_dg(d + "V", n, speed, problem.smesh);
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

} // namespace Nares
} // namespace Nextsim

#endif /* BENCHMARK_NARES_HPP */
