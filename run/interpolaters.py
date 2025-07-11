import math
import numpy as np
from scipy import interpolate
import pyproj

# Returns TOPAZ data interpolated from the data grid and coordinates to the target grid and coordinates
def topaz4_interpolate(target_lon, target_lat, data, data_x, data_y, proj_string):
    # Use pyproj and the fact that the proj4 string is embedded in the TOPAZ files
    P = pyproj.Proj(proj_string)
    target_x, target_y = P(target_lon, target_lat)

    # Mask land values and interpolate using the default griddata method (linear, as far as I can tell)
    nanmask = ~np.isnan(data.ravel())

    X, Y = np.meshgrid(data_x, data_y)
    points = np.array([X.ravel()[nanmask], Y.ravel()[nanmask]]).T
    xi = np.array([target_x.ravel(), target_y.ravel()]).T

    field = interpolate.griddata(points, data.ravel()[nanmask], xi)

    # Use griddata again to extrapolate outside the convex hull using the nearest neighbour
    nanmask = ~np.isnan(field)

    points = np.array([target_x.ravel()[nanmask], target_y.ravel()[nanmask]]).T
    xi = np.array([target_x.ravel(), target_y.ravel()]).T

    return (interpolate.griddata(points, field.ravel()[nanmask], xi, method='nearest').reshape(target_lon.shape))

# Returns bilinearly interpolated data given array of fractional indices
# 2023-03-28 Add a wrap-around for the ERA longitude. This is formally
# incorrect when this function is used for TOPAZ data, but since the target
# point would have to be out of bounds of the source, it is not so important.
def bilinear(eyes, jays, data):
    i = np.floor(eyes).astype(int)
    j = np.floor(jays).astype(int)

    fi = eyes - i
    fj = jays - j

    iwrap = (i + 1) % data.shape[1]

    return ((1 - fj) * (1 - fi) * data[j, i] +
            (1 - fj) * (fi) * data[j, iwrap] +
            (fj) * (1 - fi) * data[j + 1, i] +
            (fj) * (fi) * data[j + 1, iwrap])

# Returns ERA5 data interpolated from the data grid and coordinates to the target grid and coordinates
def era5_interpolate(target_lons, target_lats, data, data_lons, data_lats):
    target_i = (target_lons - data_lons[0]) / (data_lons[1] - data_lons[0])
    # Make sure that the index is in the range of the size of the longitude array
    target_i %= len(data_lons)

    # Latitudes are on a Gaussian grid, so we need to search a bit.
    target_j = (target_lats - data_lats[0]) / (data_lats[1] - data_lats[0])

    return bilinear(target_i, target_j, data)

# Returns the rotation angle at a given position to transform vectors from
# geographic pole coordinates to the Greenland displaced pole coordinate system
def heading_to_greenland(lat, lon):
    # The pole lies at 75˚ N, 40˚ W = -40˚ E = 320˚ E
    pole_lat = math.radians(75.0)
    pole_lon = 320.0

    delta_lon = np.radians(pole_lon - lon)

    rlat = np.radians(lat)

    # The rotation of the basis vector can be computed as the azimuth of the
    # great circle path from the location to the location of the new pole.
    return np.arctan2(math.cos(pole_lat) * np.sin(delta_lon), np.cos(rlat) * math.sin(pole_lat) - np.sin(rlat) * math.cos(pole_lat) * np.cos(delta_lon))

# Rotates the u and v velocity components by an angle given in radians
def rotate_velocities(u, v, angle):
    return (u * np.cos(angle) + v * np.sin(angle), -u * np.sin(angle) + v * np.cos(angle))
