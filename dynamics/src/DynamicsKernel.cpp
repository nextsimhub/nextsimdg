/*!
 * @file DynamicsKernel.cpp
 *
 * @date 27 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Tim Spain <piotr.minakowski@ovgu.de>
 */

#include "include/DynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection>
const std::unordered_map<std::string, ModelArray::Type> DynamicsKernel<CGdegree, DGadvection>::fieldType = {

};


} /* namespace Nextsim */
