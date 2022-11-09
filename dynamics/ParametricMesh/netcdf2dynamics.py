#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert mesh from netcdf (neXtSIM exchange grid files) file 
into nexstim dynamics format

Meshfiles can be found under the link
https://ige-meom-opendap.univ-grenoble-alpes.fr/thredds/catalog/meomopendap/extract/SASIP/grids/NH_PS/catalog.html


Given file has following numbering of the corners: 
0 bottom-left, 1 top-left, 2 top-right, 3 bottom-right 

We assume that in the 

with Nx x Ny elements

@author: Piotr Minakowski
"""

import numpy as np
import netCDF4


def check_corner_points(lon_corners,lat_corners, nx, ny):
  #check if points are consequently ordered in the middle of the domain
  #matches corresponding points of neighbourning cells
  pairs_x = [[3,0],[2,1]]
  pairs_y = [[0,1],[3,2]]
  for iy in range(ny-2):
    for ix in range(nx-2):
      for px in pairs_x:
        if not ( (lon_corners[ix  ,iy+1][px[0]] == lon_corners[ix+1,iy+1][px[1]]) 
             and (lat_corners[ix  ,iy+1][px[0]] == lat_corners[ix+1,iy+1][px[1]])  
             and (lon_corners[ix+1,iy+1][px[0]] == lon_corners[ix+2,iy+1][px[1]]) 
             and (lat_corners[ix+1,iy+1][px[0]] == lat_corners[ix+2,iy+1][px[1]])):
          print("Error with maching points")
          return False
      for py in pairs_y:
        if not ( (lon_corners[ix+1,iy+2][py[0]] == lon_corners[ix+1,iy+1][py[1]]) 
             and (lat_corners[ix+1,iy+2][py[0]] == lat_corners[ix+1,iy+1][py[1]])  
             and (lon_corners[ix+1,iy+1][py[0]] == lon_corners[ix+1,iy  ][py[1]]) 
             and (lat_corners[ix+1,iy+1][py[0]] == lat_corners[ix+1,iy  ][py[1]])):
          print("Error with maching points")
          return False
  return True




if __name__ == '__main__':


  input_file ="25km_NH.nc"


  output_file = input_file.replace(".nc",'.smesh')

  #load file
  nc = netCDF4.Dataset(input_file)


  #prepare file for save
  f = open(f"{output_file}", "w")
  

  #centers of elemenents
  lat_center = np.array(nc['plat'])
  lon_center = np.array(nc['plon'])

  #get number of elements in nx and ny
  assert(lat_center.shape == lon_center.shape)
  nx, ny = lat_center.shape

  #four corners of the grid
  lat_corners = np.array(nc['lat_corners'])
  lon_corners = np.array(nc['lon_corners'])


  assert( check_corner_points(lon_corners,lat_corners, nx, ny) )

  #Start of the mesh file
  f.write("ParametricMesh 1.0\n")
  #number of elements in x and y directions
  f.write('{0}\t{1}\n'.format(nx,ny))
  #Saving points
  for iy in range(ny):
    for ix in range(nx):
      f.write('{0}\t{1}\n'.format(lon_corners[ix,iy,0],lat_corners[ix,iy,0]))
    #add last point on the right
    f.write('{0}\t{1}\n'.format(lon_corners[ix,iy,3],lat_corners[ix,iy,3]))
  #add last row on top
  for ix in range(nx):
    f.write('{0}\t{1}\n'.format(lon_corners[ix,iy,1],lat_corners[ix,iy,1]))
  #add last point on the right in the top row
  assert(ix==nx-1)
  assert(iy==ny-1)
  f.write('{0}\t{1}\n'.format(lon_corners[ix,iy,3],lat_corners[ix,iy,3]))
  f.flush()

  #Dirichlet boundaries for land mask
  mask = np.array(nc['mask'])
  for iy in range(ny):
    for ix in range(nx):
      if mask[ix,iy] == 0: # check if we are on Ice 
        #get number of element
        no_element = ix + (nx+1)*iy
        #save corresponding edge 
        for shift in [[1,0,1],[-1,0,3],[0,1,2],[0,-1,0]]:
          #try to save boundary based on neighbour
          try:
            assert(ix+shift[0]>=0)
            assert(iy+shift[1]>=0)
            if mask[ix+shift[0],iy+shift[1]] != 0:
              f.write('{}\t{}\n'.format(no_element, shift[2]) )
          
          #except: we are on the boundary of the domain
          except:
            f.write('{}\t{}\n'.format(no_element, shift[2]) )
  
  f.close()