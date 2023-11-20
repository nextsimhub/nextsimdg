import matplotlib.pyplot as plt
import netCDF4
import numpy as np

file = 'diagnostic.nc'

# Load the basic variables
root = netCDF4.Dataset(file, "r", format="NETCDF4")
hice = np.squeeze(np.array(root.groups["data"].variables["hice"][:].data))
hsnow = np.squeeze(np.array(root.groups["data"].variables["hsnow"][:].data))
tice = np.array(root.groups["data"].variables["tice"][:].data)

# Calculate ice draught for a nicer visualisation
rho = 917
rhoSnow = 330
rhoOcean = 1025
iceDraught = (hice * rho + hsnow * rhoSnow) / rhoOcean

# Some simple diagnostics
print('hice  max: {0:0.2f}, min: {1:0.2f}, mean: {2:0.2f}'.format(hice.max(), hice.min(), hice.mean()))
print('hsnow max: {0:0.2f}, min: {1:0.2f}, mean: {2:0.2f}'.format(hsnow.max(), hsnow.min(), hsnow.mean()))

# Figure showing temperature evolution, cf,. Winton (2000) figure 2
plt.figure(1)
plt.plot([0, len(hice)], [0, 0], 'k--')
plt.plot(np.squeeze(tice[:, 0]), 'k', label="Surface")
plt.plot(np.squeeze(tice[:, 1]), label="T1")
plt.plot(np.squeeze(tice[:, 2]), label="T2")
plt.xlabel("Day of year")
plt.ylabel("Temperature [°C]")
plt.legend()
plt.show(block=False)

# Figure showing thickness evolution, cf. Winton (2000) figure 2
plt.figure(2)
plt.plot([0, len(hice)], [0, 0], 'k--')
plt.plot(hice - iceDraught, 'b', label="Ice")
plt.plot(hice + hsnow - iceDraught, 'k', label="Snow")
plt.plot(-iceDraught, 'b')
plt.xlabel("Day of year")
plt.ylabel("Height over sea level [m]")
plt.legend()
plt.show()

# TODO: Reproduce exactly Winton's figures 2 and 3