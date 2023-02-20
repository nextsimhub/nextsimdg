/*!
 * @file DummyDynamics.cpp
 *
 * @date 16 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DummyDynamics.hpp"

namespace Nextsim {

DummyDynamics::DummyDynamics()
    : hice(ModelArray::Type::DG)
    , cice(ModelArray::Type::DG)
    , hsnow(ModelArray::Type::DG)
    , uice(ModelArray::Type::CG)
    , vice(ModelArray::Type::CG)
{
}

void DummyDynamics::setData(const ModelState::DataMap& ms)
{
    hice.resize();
    cice.resize();
    hsnow.resize();
    uice.resize();
    vice.resize();


}

}
