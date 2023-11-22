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
import matplotlib.pyplot as plt


if __name__ == '__main__':


  input_file ="25km_NH.nc"
  restart_file = "init_25km_NH.nc" # get the land mask from the restart file

  output_file = input_file.replace(".nc",'_newmask.smesh')

  #load file
  nc = netCDF4.Dataset(input_file)
  restart = netCDF4.Dataset(restart_file)

  #prepare file for save
  f = open(f"{output_file}", "w")

  #Dirichlet boundaries for land mask
  dirichlet_list = [] #first collecting all dirichlet boundaries

  #We only care about mask
  #mask = np.array(nc['mask']).T
  mask = np.array(restart["data/mask"]).T

  #print(mask.shape)
  nx, ny = mask.shape

  #Start of the mesh file
  f.write("ParametricMesh 2.0\n")
  #number of elements in x and y directions
  f.write('{0}\t{1}\n'.format(nx,ny))

  x_scale = 25000
  y_scale = 25000

  #simple write mesh elements 1x1
  for iy in range(ny+1):
   for ix in range(nx+1):
      f.write('{0}\t{1}\n'.format(ix * x_scale, iy * y_scale))

  for iy in range(ny):
    for ix in range(nx):
      if mask[ix,iy] == 1: # check if we are on Ice
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
  ## write landmask
  f.write("landmask\t{}\n".format(nx*ny));
  for iy in range(ny):
    for ix in range(nx):
      f.write("{}\n".format(int(mask[ix,iy])))

  #plt.imshow(mask);plt.show() #quick check

  #write dirichlet boundaries into file
  f.write("dirichlet\t{}\n".format(len(dirichlet_list)))
  for dirichlet in dirichlet_list:
    f.write('{}\t{}\n'.format(dirichlet[0], dirichlet[1]))

  f.write("periodic 0")
  f.close()
