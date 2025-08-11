import math

import numpy as np
import pyproj
from scipy import interpolate


def topaz4_interpolate(target_lon, target_lat, data, data_x, data_y, proj_string):
    """
    Returns TOPAZ data interpolated from the data grid and coordinates to the target grid and coordinates.

    :param target_lon: Longitudes of the target grid
    :param target_lat: Latitudes of the target grid
    :param data: Data to interpolate
    :param data_x: X coordinates of the data grid
    :param data_y: Y coordinates of the data grid
    :param proj_string: Map projection (proj4 string) of the data grid
    :return: Interpolated data on the target grid
    """
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

    return (interpolate.griddata(points, field.ravel()[nanmask], xi, method="nearest").reshape(target_lon.shape))

def bilinear(eyes, jays, data):
    """"
    Returns bilinearly interpolated data given an array of fractional indices. Supports periodic boundary conditions.

    :param eyes: Fractional indices along the x-axis
    :param jays: Fractional indices along the y-axis
    :param data: Data to interpolate
    """
    i = np.floor(eyes).astype(int)
    j = np.floor(jays).astype(int)

    fi = eyes - i
    fj = jays - j

    iwrap = (i + 1) % data.shape[1]

    return ((1 - fj) * (1 - fi) * data[j, i] +
            (1 - fj) * (fi) * data[j, iwrap] +
            (fj) * (1 - fi) * data[j + 1, i] +
            (fj) * (fi) * data[j + 1, iwrap])

def era5_interpolate(target_lons, target_lats, data, data_lons, data_lats):
    """
    Returns ERA5 data interpolated from the data grid and coordinates to the target grid and coordinates.

    :param target_lons: Longitudes of the target grid
    :param target_lats:  Latitudes of the target grid
    :param data: Data to interpolate
    :param data_lons: Longitudes of the data grid
    :param data_lats: Latitudes of the data grid
    :return: Interpolated data on the target grid
    """
    target_i = (target_lons - data_lons[0]) / (data_lons[1] - data_lons[0])
    # Make sure that the index is in the range of the size of the longitude array
    target_i %= len(data_lons)

    # Latitudes are on a Gaussian grid, so we need to search a bit.
    target_j = (target_lats - data_lats[0]) / (data_lats[1] - data_lats[0])

    return bilinear(target_i, target_j, data)

def heading_to_greenland(lat, lon):
    """
    Returns the rotation angle at a given position to transform vectors from geographic pole coordinates to the Greenland displaced pole coordinate system.

    :param lat: Latitude of the location
    :param lon: Longitude of the location
    :return: The rotation of the basis vector, computed as the azimuth of the great circle path from the location to the location of the new pole.
    """
    # The pole lies at 75˚ N, 40˚ W = -40˚ E = 320˚ E
    pole_lat = math.radians(75.0)
    pole_lon = 320.0

    delta_lon = np.radians(pole_lon - lon)

    rlat = np.radians(lat)

    # The rotation of the basis vector can be computed as the azimuth of the
    # great circle path from the location to the location of the new pole.
    return np.arctan2(math.cos(pole_lat) * np.sin(delta_lon), np.cos(rlat) * math.sin(pole_lat) - np.sin(rlat) * math.cos(pole_lat) * np.cos(delta_lon))

def rotate_velocities(u, v, angle):
    """
    Rotates the u and v velocity components by an angle given in radians.

    :param u: U-component of the velocity vector
    :param v: V-component of the velocity vector
    :param angle: The angle to rotate by, in radians
    :return: A tuple of the rotated u and v components.
    """
    return (u * np.cos(angle) + v * np.sin(angle), -u * np.sin(angle) + v * np.cos(angle))

def rotate_pole_to_greenland(lat, lon):
    """
    Rotates the mesh such that the singularities are in Greenland / Antarctica at 75N / 40W and 75S / 140E
    This is a copy of ParametricMesh::RotatePoleToGreenland in dynamics/src/include/ParametricMesh.hpp. It's here for testing purposes only.

    :param lat: Latitudes of the mesh
    :param lon: Longitudes of the mesh
    :return: A tuple of the rotated latitudes and longitudes.
    """
    latr = np.deg2rad(lat)
    lonr = np.deg2rad(lon)

    x = np.cos(latr) * np.cos(lonr)
    y = np.cos(latr) * np.sin(lonr)
    z = np.sin(latr)

    aw = np.deg2rad(40.0)
    x1 = np.cos(aw) * x - np.sin(aw) * y
    y1 = np.sin(aw) * x + np.cos(aw) * y
    z1 = z

    bw = np.deg2rad(15.0)
    x2 = np.cos(bw) * x1 - np.sin(bw) * z1
    y2 = y1
    z2 = np.sin(bw) * x1 + np.cos(bw) * z1

    return (np.rad2deg(np.arcsin(z2)), np.rad2deg(np.arctan2(y2, x2)))
