/*!
 * @file VectorManupulations.hpp
 * @date 9 Juli 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __VECTORMANIPULATIONS_HPP
#define __VECTORMANIPULATIONS_HPP

/*!
 * This namespace combines all actions on vectors
 * that only depend on the mesh but not on further
 * knowledge, such as the discretization type or
 * the equation to be solved.
 *
 * Examples are the handling of the boundaries
 */

#include "ParametricMesh.hpp"
#include "cgVector.hpp"

namespace Nextsim {

namespace VectorManipulations {

    template <int CG> void CGAveragePeriodic(const ParametricMesh& smesh, CGVector<CG>& v)
    {
        // the two segments bottom, right, top, left, are each processed in parallel
        for (size_t seg = 0; seg < smesh.periodic.size(); ++seg) {
            // #pragma omp parallel for
            for (size_t i = 0; i < smesh.periodic[seg].size(); ++i) {

                const size_t ptype = smesh.periodic[seg][i][0];
                const size_t eid_lb = smesh.periodic[seg][i][2];
                const size_t eid_rt = smesh.periodic[seg][i][1];

                size_t ix_lb = eid_lb % smesh.nx;
                size_t iy_lb = eid_lb / smesh.nx;
                size_t i0_lb = (CG * smesh.nx + 1) * CG * iy_lb
                    + CG * ix_lb; // lower/left index in left/bottom element
                size_t ix_rt = eid_rt % smesh.nx;
                size_t iy_rt = eid_rt / smesh.nx;
                size_t i0_rt = (CG * smesh.nx + 1) * CG * iy_rt
                    + CG * ix_rt; // lower/left index in right/top element

                if (ptype == 0) // X-edge, bottom/top
                {
                    for (size_t j = 0; j <= CG; ++j) {
                        v(i0_lb + j)
                            = 0.5 * (v(i0_lb + j) + v(i0_rt + CG * (CG * smesh.nx + 1) + j));
                        v(i0_rt + CG * (CG * smesh.nx + 1) + j) = v(i0_lb + j);
                    }
                } else if (ptype == 1) // Y-edge, left/right
                {
                    for (size_t j = 0; j <= CG; ++j) {
                        const size_t i1 = i0_lb + j * (CG * smesh.nx + 1);
                        const size_t i2 = i0_rt + CG + j * (CG * smesh.nx + 1);
                        v(i1) = 0.5 * (v(i1) + v(i2));
                        v(i2) = v(i1);
                    }
                } else
                    abort();
            }
        }
    }

    template <int CG> void CGAddPeriodic(const ParametricMesh& smesh, CGVector<CG>& v)
    {
        // the two segments bottom, right, top, left, are each processed in parallel
        for (size_t seg = 0; seg < smesh.periodic.size(); ++seg) {
            // #pragma omp parallel for
            for (size_t i = 0; i < smesh.periodic[seg].size(); ++i) {

                const size_t ptype = smesh.periodic[seg][i][0];
                const size_t eid_lb = smesh.periodic[seg][i][2];
                const size_t eid_rt = smesh.periodic[seg][i][1];

                size_t ix_lb = eid_lb % smesh.nx;
                size_t iy_lb = eid_lb / smesh.nx;
                size_t i0_lb = (CG * smesh.nx + 1) * CG * iy_lb
                    + CG * ix_lb; // lower/left index in left/bottom element
                size_t ix_rt = eid_rt % smesh.nx;
                size_t iy_rt = eid_rt / smesh.nx;
                size_t i0_rt = (CG * smesh.nx + 1) * CG * iy_rt
                    + CG * ix_rt; // lower/left index in right/top element

                //      problem: wenn man bis CG geht, dann werden die ecken je doppelt addiert...
                //               wenn man eins vorher aufhoert, dann ist in der ecke was falsch?
                if (ptype == 0) // X-edge, bottom/top
                {
                    for (size_t j = 0; j < CG; ++j) {
                        v(i0_lb + j) += v(i0_rt + CG * (CG * smesh.nx + 1) + j);
                        v(i0_rt + CG * (CG * smesh.nx + 1) + j) = v(i0_lb + j);
                    }
                } else if (ptype == 1) // Y-edge, left/right
                {
                    for (size_t j = 0; j < CG; ++j) {
                        const size_t i1 = i0_lb + j * (CG * smesh.nx + 1);
                        const size_t i2 = i0_rt + CG + j * (CG * smesh.nx + 1);

                        v(i1) += v(i2);
                        v(i2) = v(i1);
                    }
                } else
                    abort();
            }
        }
    }

}

}

#endif // __VECTORMANIPULATIONS_HPP
