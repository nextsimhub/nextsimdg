# Standalone dynamics benchmarks

Executables that run a dynamics kernel of nextSIM-DG without the rest of the model, on the
sea-ice benchmark of Mehlmann et al. (2021): a 512 km × 512 km box (with one layer of land cells
around it), ice concentration 1, ice height `0.3 + 0.005 (sin(6e-5 x) + sin(3e-5 y))`, a
circular ocean current and a cyclone moving to the north east. The forcing is the same as in
the modules `BenchmarkOcean` and `BenchmarkAtmosphere` of the full model.

| Executable                    | Kernel                  | Rheology / momentum                     |
|-------------------------------|-------------------------|-----------------------------------------|
| `benchmark_mehlmann_mevp`     | `MEVPDynamicsKernel`    | mEVP, CG velocity                       |
| `benchmark_mehlmann_bbm`      | `BBMDynamicsKernel`     | brittle Bingham-Maxwell, CG velocity    |
| `benchmark_mehlmann_ldg_mevp` | `LDGMEVPDynamicsKernel` | mEVP, LDG (velocity space `LDGVelocity`) |

Build with the CMake option `BUILD_BENCHMARKS`:

```
cmake -B build -S . -DDynamicsType=DQ1 -DLDGVelocity=DQ1 -DBUILD_BENCHMARKS=ON
make -C build benchmark_mehlmann_mevp benchmark_mehlmann_bbm benchmark_mehlmann_ldg_mevp
```

The advection space is that of `DynamicsType`, the CG velocity degree and the LDG velocity
space follow `DynamicsType` and `LDGVelocity`, as for the model.

Options (defaults in brackets):

```
--nx N          cells in each direction (64)
--dt seconds    time step (120)
--days d        simulated time (2)
--output hours  interval of the VTK output, 0: no output (6)
--nsteps n      mEVP / BBM sub-iterations per time step (120)
--dir name      output directory (name of the executable)
```

The output directory contains VTK files of the ice height `H`, concentration `A`, damage `D`
(BBM only) as DG fields and of the cell averages of the velocity `u`, `v` and of the shear.
