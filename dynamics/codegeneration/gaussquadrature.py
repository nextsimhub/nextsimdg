## module used in the code-generation scripts

import numpy as np


### Gauss quadrature
gausspoints = np.array([
    [0.5,0,0,0],
    [0.5-np.sqrt(1./12.),0.5+np.sqrt(1./12.),0,0],
    [0.5-np.sqrt(3./20.),0.5,0.5+np.sqrt(3./20.),0],
    [0.5-0.5*np.sqrt(3./7.+2./7.*np.sqrt(6./5.)),
     0.5-0.5*np.sqrt(3./7.-2./7.*np.sqrt(6./5.)),
     0.5+0.5*np.sqrt(3./7.-2./7.*np.sqrt(6./5.)),
     0.5+0.5*np.sqrt(3./7.+2./7.*np.sqrt(6./5.))]])

gaussweights = np.array([
    [1.0,0,0,0],
    [0.5,0.5,0,0],
    [5./18.,8./18.,5./18.,0],
    [(18.-np.sqrt(30.0))/72.,
     (18.+np.sqrt(30.0))/72.,
     (18.+np.sqrt(30.0))/72.,
     (18.-np.sqrt(30.0))/72.]])

### LagrangeQuadrature (Trapez / Simpson)
lagrangepoints = np.array([
    [0.0,1.0,0,0],
    [0.0,0.5,1,0]])

lagrangeweights = np.array([
    [0.5,0.5,0.0],
    [1./6.,2./3.,1./5.]])

def sanitycheck_gauss():
    for p in range(10):    # integrate x^p for p=0,1,2,...,9
        ex = 1.0/(1.0+p)
#        print("Integrate x^{0} over [0,1]: ".format(p),end='')

        for q in range(4): # gauss 1,2,3,4
            gi = 0.0
            for k in range(q+1):
                gi = gi + gaussweights[q,k] * gausspoints[q,k]**p
            if np.fabs(gi-ex)>1.e-14:
                assert p>=2*(q+1), 'Gauss rule should be exact'
