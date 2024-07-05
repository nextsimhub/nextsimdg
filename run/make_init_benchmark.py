from make_init_base import initMaker
from math import sin

# Creates initial conditions for the Mehlmann et al. (2021) benchmark case, at 2, 4, 8, and 16 km resolutions.

config_contents = """
[model]
init_file = init_benchmark_{nfirst}x{nsecond}.nc
start = 2023-01-01T00:00:00Z
stop = 2023-01-01T01:00:00Z
time_step = P0-0T00:02:00

[Modules]
DiagnosticOutputModule = Nextsim::ConfigOutput
DynamicsModule = Nextsim::MEVPDynamics
IceThermodynamicsModule = Nextsim::DummyIceThermodynamics
LateralIceSpreadModule = Nextsim::DummyIceSpread
AtmosphereBoundaryModule = Nextsim::BenchmarkAtmosphere
OceanBoundaryModule = Nextsim::BenchmarkOcean

[ConfigOutput]
period = P0-0T04:00:00
start = 2010-01-01T12:00:00Z
field_names = hsnow,hice,cice,u,v
filename = benchmark_{nfirst}x{nsecond}_gpu.diagnostic.nc
"""

# Domain size [km]
L = 512
for res in [0.125, 0.25, 0.5, 1, 2, 4]:

    nfirst = int(L / res)
    nsecond = int(L / res)
    nLayers = 1

    with open(f"config_benchmark_{nfirst}x{nsecond}.cfg", 'w') as f:
        f.write(config_contents.format(nfirst=nfirst, nsecond=nsecond))

    fname = f"init_benchmark_{nfirst}x{nsecond}.nc"
    print("Producing file", fname)

    initializer = initMaker(fname, nfirst, nsecond, nLayers, res*1e3)
    # The model expects everything in metres, while the benchmark problem in Mehlman et al. (2021) is defined in km.

    # Ice everywhere and all boundaries closed
    initializer.mask[:, :] = 1.
    initializer.mask[0, :] = 0.
    initializer.mask[-1, :] = 0.
    initializer.mask[:, 0] = 0.
    initializer.mask[:, -1] = 0.

    # Uniform concentration of 100%
    initializer.cice[:, :] = 1.

    # Loop over ice thickness to construct the initial conditions. This should be a pattern of undulating ice.
    for ix in range(nfirst):
        x = ix * res
        for iy in range(nsecond):
            y = iy * res
            initializer.hice[ix, iy] = 0.3 + 0.005 * (sin(60e-3 * x) + sin(30e-3 * y))

    # All other variables are zero or not needed

    # The file is written when initializer goes out of scope
