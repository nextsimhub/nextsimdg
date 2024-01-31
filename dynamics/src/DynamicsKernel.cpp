/*!
 * @file DynamicsKernel.cpp
 *
 * @date Jan 31, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DynamicsKernel.hpp"

namespace Nextsim {

template class DynamicsKernel<1, 1>;
template class DynamicsKernel<1, 2>;
template class DynamicsKernel<1, 3>;

template class DynamicsKernel<2, 1>;
template class DynamicsKernel<2, 2>;
template class DynamicsKernel<2, 3>;
}


