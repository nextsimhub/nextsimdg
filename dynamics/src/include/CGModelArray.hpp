/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CGMODELARRAY_HPP
#define CGMODELARRAY_HPP

#include "Interpolations.hpp"
#include "cgVector.hpp"

#include "include/ModelArray.hpp"

namespace Nextsim {

class CGModelArray {
public:
    template <int CG> static CGVector<CG>& ma2cg(const ModelArray& ma, CGVector<CG>& cg)
    {
        if (ma.getType() != ModelArray::Type::CG) {
            /*
             * Create a ParametricMesh with the correct x and y dimensions, the
             * only members used by the Interpolations functions. Constructed
             * with Cartesian coordinates, but the coordinate system is not used.
             */
            ParametricMesh smesh(Nextsim::CARTESIAN); // The coordinate system is unimportant here
            smesh.nx = ma.dimensions()[0];
            smesh.ny = ma.dimensions()[1];
            // Assume the data is compatible with a DG0 array
            DGVector<1> asDG(smesh);
            asDG = ma.data().matrix();
            Interpolations::DG2CG(smesh, cg, asDG);
        } else {
            // CG to CG. Assume the CG degrees are equal
            cg = ma.data().matrix();
        }
        return cg;
    }

    template <int CG> static ModelArray& cg2ma(const CGVector<CG>& cg, ModelArray& ma)
    {
        ma.setData(cg.data());
        return ma;
    }

    template <int CG> static ModelArray::MultiDim cgDimensions(const ModelArray::MultiDim& hDims)
    {
        ModelArray::MultiDim cgDims(hDims);
        for (size_t i = 0; i < cgDims.size(); ++i) {
            cgDims[i] *= CG;
            ++cgDims[i];
        }
        return cgDims;
    }
};
}

#endif /* CGMODELARRAY_HPP */
