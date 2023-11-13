/*!
 * @file DynamicsKernel.cpp
 *
 * @date 27 Mar 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Piotr Minakowski <piotr.minakowski@ovgu.de>
 */

#include "include/MEVPDynamicsKernel.hpp"

namespace Nextsim {

template <int CGdegree, int DGadvection>
const std::unordered_map<std::string, ModelArray::Type> MEVPDynamicsKernel<CGdegree, DGadvection>::fieldType = {

};


} /* namespace Nextsim */
