import netCDF4
import glob
import numpy as np
import matplotlib.pyplot as plt

first = True
for file in sorted(glob.glob('diagnostic.*.nc')):

    root = netCDF4.Dataset(file, "r", format="NETCDF4")

    time_var = root.groups['metadata'].groups["time"].variables["time"]
    np_time = np.datetime64(netCDF4.num2date(time_var[:], time_var.units).isoformat())

    if first:
        time = np.array(np_time)
        hice = np.array(root.groups["data"].variables["hice"][0].data)
        hsnow = np.array(root.groups["data"].variables["hsnow"][0].data)
        tice = np.array(root.groups["data"].variables["tice"][0].data)
        first = False
    else:
        time = np.append(time, np_time)
        hice = np.append(hice, root.groups["data"].variables["hice"][0].data)
        hsnow = np.append(hsnow, root.groups["data"].variables["hsnow"][0].data)
        tice = np.vstack((tice, root.groups["data"].variables["tice"][0].data))

end = len(time) - 1

plt.figure(1)
plt.plot([time[0], time[end]], [0, 0])
plt.plot(time, tice[:, 0])
# plt.plot(time, tice[:, 1])
# plt.plot(time, tice[:, 2])
plt.show(block=False)

plt.figure(2)
plt.plot(time, hice)
plt.plot(time, hice+hsnow)
ax = plt.gca()
ax.set_ylim([0, None])
plt.show()
