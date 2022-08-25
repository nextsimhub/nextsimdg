#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 14:53:04 2022

Creates a simple mesh of size [0,Lx] x [0,Ly]
with Nx x Ny elements

The mesh is slightly distorted

@author: richter
"""

# Size of the domain
Lx =  512000
Ly =  512000

# Number of elements
nx =  32
ny =  32

nx = ny = 512

import numpy as np

def distx(x,y):
    return np.sin(3.0*np.pi*x/Lx)*np.sin(2.0*np.pi*y/Ly)

def disty(x,y):
    return np.sin(np.pi*x/Lx)*np.sin(3.0*np.pi*y/Ly)



sx = 0.05*Lx
sy = 0.05*Ly

sx,sy = 0,0

f = open(f"distortedrectangle_{nx}x{ny}.smesh", "w")
f.write("SasipMesh 1.0\n")
f.write('{0}\t{1}\n'.format(nx,ny))
for y in np.linspace(0,Ly,ny+1):
    for x in np.linspace(0,Lx,nx+1):
        f.write('{0}\t{1}\n'.format(x+sx*distx(x,y),y+sy*disty(x,y)))

f.close()
