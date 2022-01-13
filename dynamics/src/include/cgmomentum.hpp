/*----------------------------   cgmomentum.hpp     ---------------------------*/
#ifndef __cgmomentum_HPP
#define __cgmomentum_HPP
/*----------------------------   cgmomentum.hpp     ---------------------------*/

#include "cgvector.hpp"
#include "dgvector.hpp"

namespace Nextsim {

class CGMomentum {
private:
public:
    /*!
     * The following functions take care of the interpolation and projection
     * between CG and DG functions
     */
    //! Projects a CG function to a DG function
    template <int CG, int DG>
    void ProjectCGToDG(const Mesh& mesh, CellVector<DG>& dg, const CGVector<CG>& cg);

    //! Projects the symmetric gradient of the CG2 velocity into the DG1 space
    template <int CG, int DG>
    void ProjectCG2VelocityToDG1Strain(const Mesh& mesh,
        CellVector<DG>& E11, CellVector<DG>& E12, CellVector<DG>& E22,
        const CGVector<CG>& vx, const CGVector<CG>& vy);

    /*!
   * Evaluates (S, nabla phi) and adds it to tx/ty - Vector
   */
    template <int CG, int DG>
    void AddStressTensor(const Mesh& mesh, double scale,
        CGVector<CG>& tx, CGVector<CG>& ty,
        const CellVector<DG>& S11, const CellVector<DG>& S12, const CellVector<DG>& S22) const;

    //! Sets the vector to zero along the boundary
    template <int CG>
    void DirichletZero(const Mesh& mesh, CGVector<CG>& v) const;
};

}

/*----------------------------   cgmomentum.hpp     ---------------------------*/
/* end of #ifndef __cgmomentum_HPP */
#endif
/*----------------------------   cgmomentum.hpp     ---------------------------*/
