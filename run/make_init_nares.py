from make_init_base import initMaker
from math import sin
import numpy as np

# Creates initial conditions for the Mehlmann et al. (2021) nares case, at 2, 4, 8, and 16 km resolutions.

'''
   *  (0,120)-+- (80,120)                       / (240,120)   H2 = 120
   *          |     |    \                 /   +     |
   *          |     |    (130,80) -(170,80)    |     |        H1 = 80
   *  + ------------+---- +-------------+------+-----|
   *          |     |    (130,40) -(170,40)    |     |        H0 = 40
   *          |     |   /                  \   +     |
   *  (0,0)---+-- (80,0)                        \  (240,0)   0
   *
   * 1st: y < 40/50*(x-80)      <=> 5/4 y + 80 < x
   * 2nd: y < 40-40/70*(x-170)  <=> (40-y)*7/4+170 > x
'''

# Domain size [km]
Ly = 120
Lx = 2*Ly

for N in [32, 64, 128]:  # number of elements in y-direction

    ny = N
    nx = 2*N
    hx = Lx/nx
    hy = Ly/ny

    assert hx==hy
    
    nLayers = 1

    fname = f"init_nares_{nx}x{ny}.nc"

    initializer = initMaker(fname, ny, nx, nLayers, hx*1e3, checkZeros=False)
    # The model expects everything in metres, while the nares problem in Mehlman et al. (2021) is defined in km.

    # Ice everywhere and all boundaries closed
    initializer.mask[:, :] = 1.
    initializer.mask[0, :] = 0.
    initializer.mask[-1, :] = 0.
#    initializer.mask[:, 0] = 0.  # left/ right open
#    initializer.mask[:, -1] = 0. # 

    # narrowing
    for iy in range(ny):
        for ix in range(nx):
            x = (0.5+ix)*hx
            y = (0.5+iy)*hy

            dy = min(y,1.*Ly-y)
            if dy < 20 and x > 80+5./4. * dy and x < 170. + 7./4.*(40.-dy):
                initializer.mask[iy,ix] = 0
                    
                
    
    # Uniform concentration of 100%
    initializer.cice[:, :] = 1.

    # Loop over ice thickness to construct the initial conditions. This should be a pattern of undulating ice.
    for ix in range(nx):
        x = (0.5+ix) * hx * 1.e3
        for iy in range(ny):
            y = (0.5 + iy) * hy * 1.e3
            initializer.hice[iy, ix] = 0.3 + 0.005 * (sin(60e-3 * x) + sin(30e-3 * y))

    initializer.damage[:, :] = 1.

    # All other variables are zero or not needed

    # The file is written when initializer goes out of scope
