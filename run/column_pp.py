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
print('tsurf max: {0:0.2f}, min: {1:0.2f}, mean: {2:0.2f}'.format(tice[:, 0].max(), tice[:, 0].min(), tice[:, 0].mean()))
print('t1    max: {0:0.2f}, min: {1:0.2f}, mean: {2:0.2f}'.format(tice[:, 1].max(), tice[:, 1].min(), tice[:, 1].mean()))
print('t2    max: {0:0.2f}, min: {1:0.2f}, mean: {2:0.2f}'.format(tice[:, 2].max(), tice[:, 2].min(), tice[:, 2].mean()))

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
plt.show(block=False)

# TODO: Reproduce exactly Winton's figures 2 and 3

fig, (ax1, ax2, ax3) = plt.subplots(3)

ax1.plot(np.linspace(0, 13, 365), hice, 'k')
ax1.plot(np.linspace(0, 13, 365), hice + hsnow, 'k')
ax1.set_xlim([0, 13])
ax1.set_ylim([2.4, 3.7])
ax1.grid(True)

ax2.plot(np.linspace(0, 13, 365), np.squeeze(tice[:, 1]))
ax2.set_xlim([0, 13])
ax2.set_ylim([-15, 0])
ax2.grid(True)

ax3.plot([0, len(hice)], [-1.8, -1.8], 'k--')
ax3.plot(np.linspace(0, 13, 365), np.squeeze(tice[:, 2]))
ax3.set_xlim([0, 13])
ax3.set_ylim([-6, -1])
ax3.grid(True)

plt.show()
