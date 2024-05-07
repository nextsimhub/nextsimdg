from make_init_base import initMaker

nfirst = 1
nsecond = 1
nLayers = 3
resolution = 50

initializer = initMaker("init_column.nc", nfirst, nsecond, nLayers, resolution)

ice_salinity = 5  # should match Ice::s in constants.hpp
mu: float = -0.055  # should match Water::mu in constants.hpp
ocean_temperature = -1.54
ocean_salinity = ocean_temperature / mu

initializer.mask[:, :] = 1
initializer.cice[:, :] = 1.
initializer.hice[:, :] = 2.
initializer.hsnow[:, :] = 0.3
initializer.sss[:, :] = ocean_salinity
initializer.sst[:, :] = ocean_temperature
initializer.tice[:, :, :] = ice_salinity * mu