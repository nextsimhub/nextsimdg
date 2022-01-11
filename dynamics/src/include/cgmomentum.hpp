/*----------------------------   cgmomentum.hpp     ---------------------------*/
#ifndef __cgmomentum_HPP
#define __cgmomentum_HPP
/*----------------------------   cgmomentum.hpp     ---------------------------*/

#include "cgvector.hpp"
#include "dgvector.hpp"

namespace Nextsim {

class CGMomentum {
private:
    CGVector<2> vx, vy;

public:
    /*!
     * Initializes the velocity vector according to the mesh
     */
    void ReInit(const Mesh& mesh)
    {
        vx.resize_by_mesh(mesh);
        vy.resize_by_mesh(mesh);
    }

    /*!
     * The following functions take care of the interpolation and projection
     * between CG and DG functions
     */
    //! Projects a CG function to a DG function
    template <int CG, int DG>
    void ProjectCGToDG(const Mesh& mesh, CellVector<DG>& dg, const CGVector<CG>& cg);

    //! Projects the symmetric gradient of the CG2 velocity into the DG1 space
    void ProjectCG2VelocityToDG1Strain(const Mesh& mesh,
        CellVector<1>& E11, CellVector<1>& E12, CellVector<1>& E22,
        const CGVector<2>& vx, const CGVector<2>& vy);
};

}

/*----------------------------   cgmomentum.hpp     ---------------------------*/
/* end of #ifndef __cgmomentum_HPP */
#endif
/*----------------------------   cgmomentum.hpp     ---------------------------*/
