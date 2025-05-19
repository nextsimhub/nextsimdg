/*!
 * @file dgVectorHolder.hpp
 *
 * @date Feb 4, 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef DGVECTORHOLDER_HPP
#define DGVECTORHOLDER_HPP

#include "dgVector.hpp"
#include "include/ModelArray.hpp"

namespace Nextsim {

template <int DG> class DGVectorHolder {
public:
    using EigenDGVector = typename DGVector<DG>::EigenDGVector;
    DGVectorHolder()
        : ref(nullptr)
    {
    }
    DGVectorHolder(ModelArray& ma)
        : DGVectorHolder(static_cast<ModelArray::DataType&>(ma))
    {
    }
    DGVectorHolder(ModelArray::DataType& madt)
        : DGVectorHolder(reinterpret_cast<EigenDGVector&>(madt))
    {
    }
    DGVectorHolder(EigenDGVector& edgv)
        : ref(&edgv)
    {
    }
    DGVectorHolder(DGVector<DG>& dgv)
        : ref(&dgv)
    {
    }

    operator DGVector<DG>&() { return reinterpret_cast<DGVector<DG>&>(*ref); }
    operator const DGVector<DG>&() const { return reinterpret_cast<const DGVector<DG>&>(*ref); }

    double& operator()(size_t i, size_t j) { return (*ref)(i, j); }

    void zero() { ref->setZero(); }

private:
    EigenDGVector* ref;
};
}

#endif /* DGVECTORHOLDER_HPP */
