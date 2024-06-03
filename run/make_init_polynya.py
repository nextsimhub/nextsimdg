from make_init_base import initMaker

# Creates initial conditions for the Bjornsson et al. (2001) polynya case

# Domain size [km]
x = 100
y = 50
res = 2

nfirst = int(y / res)
nsecond = int(x / res)
nLayers = 3

fname = f"init_polynya.nc"

# The model expects everything in metres
initializer = initMaker(fname, nfirst, nsecond, nLayers, res*1e3)

# Ice everywhere and all boundaries closed, except the x = 100 km end
initializer.mask[:, :] = 1.
initializer.mask[0, :] = 0.
initializer.mask[-1, :] = 0.
initializer.mask[:, 0] = 0.

# Uniform concentration of 90%
initializer.cice[:, :] = 0.9

# Uniform thickness of 20 cm
initializer.hice[:, :] = 0.2

# Ice and ocean temperature and salinity at the freezing point
ice_salinity = 5  # should match Ice::s in constants.hpp
mu: float = -0.055  # should match Water::mu in constants.hpp
ocean_temperature = -1.54
ocean_salinity = ocean_temperature / mu

initializer.sss[:, :] = ocean_salinity
initializer.sst[:, :] = ocean_temperature
initializer.tice[:, :, :] = ice_salinity * mu

# All other variables are zero or not needed

# The file is written when initializer goes out of scope
