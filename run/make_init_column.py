from make_init_base import initMaker

nfirst = 1
nsecond = 1
nLayers = 3
resolution = 50

initializer = initMaker("init_column.nc", nfirst, nsecond, nLayers, resolution)

ice_salinity = 5  # should match Ice::s in constants.hpp
ocean_salinity = 32.
mu = -0.055  # should match Water::mu in constants.hpp

initializer.mask[:, :] = 1
initializer.cice[:, :] = 1.
initializer.hice[:, :] = 2.
initializer.hsnow[:, :] = 0.3
initializer.sss[:, :] = ocean_salinity
initializer.sst[:, :] = ocean_salinity * mu
initializer.tice[:, :, :] = ice_salinity * mu
