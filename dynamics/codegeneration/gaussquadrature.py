## module used in the code-generation scripts

import numpy as np

### Gauss quadrature
gausspoints = np.array(
    [
        [0.5, 0, 0, 0],
        [0.5 - np.sqrt(1.0 / 12.0), 0.5 + np.sqrt(1.0 / 12.0), 0, 0],
        [0.5 - np.sqrt(3.0 / 20.0), 0.5, 0.5 + np.sqrt(3.0 / 20.0), 0],
        [
            0.5 - 0.5 * np.sqrt(3.0 / 7.0 + 2.0 / 7.0 * np.sqrt(6.0 / 5.0)),
            0.5 - 0.5 * np.sqrt(3.0 / 7.0 - 2.0 / 7.0 * np.sqrt(6.0 / 5.0)),
            0.5 + 0.5 * np.sqrt(3.0 / 7.0 - 2.0 / 7.0 * np.sqrt(6.0 / 5.0)),
            0.5 + 0.5 * np.sqrt(3.0 / 7.0 + 2.0 / 7.0 * np.sqrt(6.0 / 5.0)),
        ],
    ]
)

gaussweights = np.array(
    [
        [1.0, 0, 0, 0],
        [0.5, 0.5, 0, 0],
        [5.0 / 18.0, 8.0 / 18.0, 5.0 / 18.0, 0],
        [
            (18.0 - np.sqrt(30.0)) / 72.0,
            (18.0 + np.sqrt(30.0)) / 72.0,
            (18.0 + np.sqrt(30.0)) / 72.0,
            (18.0 - np.sqrt(30.0)) / 72.0,
        ],
    ]
)

### LagrangeQuadrature (Trapez / Simpson)
lagrangepoints = np.array([[0.0, 1.0, 0, 0], [0.0, 0.5, 1, 0]])

lagrangeweights = np.array([[0.5, 0.5, 0.0], [1.0 / 6.0, 2.0 / 3.0, 1.0 / 5.0]])


def sanitycheck_gauss():
    for p in range(10):  # integrate x^p for p=0,1,2,...,9
        ex = 1.0 / (1.0 + p)
        #        print("Integrate x^{0} over [0,1]: ".format(p),end='')

        for q in range(4):  # gauss 1,2,3,4
            gi = 0.0
            for k in range(q + 1):
                gi = gi + gaussweights[q, k] * gausspoints[q, k] ** p
            if (np.fabs(gi - ex) > 1.0e-14) and (p < 2 * (q + 1)):
                msg = "Gauss rule should be exact"
                raise ValueError(msg)
