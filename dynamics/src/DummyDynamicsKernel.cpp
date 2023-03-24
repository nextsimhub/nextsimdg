/*!
 * @file DummyDynamicsKernel.cpp
 *
 * @date 17 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DummyDynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection>
const std::unordered_map<std::string, ModelArray::Type> DummyDynamicsKernel<CGdegree, DGadvection>::fieldType = {

};


} /* namespace Nextsim */
