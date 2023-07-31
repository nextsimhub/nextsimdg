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
@author: Tim Spain
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
  restart_file = "init_25km_NH.nc" # get the land mask from the restart file

  output_file = input_file.replace(".nc",'_newmask.smesh')

  #load file
  nc = netCDF4.Dataset(input_file)
  restart = netCDF4.Dataset(restart_file)

  #prepare file for save
  f = open(f"{output_file}", "w")
  

  #centers of elemenents
  lat_center = np.array(nc['plat'])
  lon_center = np.array(nc['plon'])

  #get number of elements in nx and ny
  assert(lat_center.shape == lon_center.shape)
  ny, nx = lat_center.shape

  #Start of the mesh file
  f.write("ParametricMesh 2.0\n")
  #number of elements in x and y directions
  f.write('{0}\t{1}\n'.format(nx,ny))
  #Saving points
  x_scale = 25000
  y_scale = 25000
  for j in range(ny+1):
      for i in range(nx+1):
          f.write(f"{i*x_scale}\t{j*y_scale}\n")
  f.flush()

  # Land-sea mask  
  mask = np.array(restart['data/mask'])
  f.write(f"landmask\t{nx*ny}\n")
  # write the mask
  for j in range(ny):
    for i in range(nx):
      f.write(f"{int(mask[j, i])}\n")
  f.flush()
  
  #Dirichlet boundaries for land mask
  dirichlet_list = [] #first collecting all dirichlet boundaries
  for iy in range(ny):
    for ix in range(nx):
      if mask[iy,ix] == 1: # check if we are on ocean 
        #get number of element
        no_element = ix + nx*iy
        #save corresponding edge 
        for shift in [[1,0,1],[-1,0,3],[0,1,2],[0,-1,0]]:
          #try to save boundary based on neighbour
          try:
            assert(ix+shift[0]>=0)
            assert(iy+shift[1]>=0)
            if mask[ix+shift[0],iy+shift[1]] != 1:
              dirichlet_list.append([no_element, shift[2]])
          #except: we are on the boundary of the domain
          except:
            dirichlet_list.append([no_element, shift[2]])

  #write dirichlet boundaries into file
  f.write("dirichlet\t{}\n".format(len(dirichlet_list)))
  for dirichlet in dirichlet_list:
    f.write('{}\t{}\n'.format(dirichlet[0], dirichlet[1]))

  f.write("periodic 0")
  f.close()