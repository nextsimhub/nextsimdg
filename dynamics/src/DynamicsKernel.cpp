/*!
 * @file DynamicsKernel.cpp
 *
 * @date Jan 31, 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/DynamicsKernel.hpp"

namespace Nextsim {

template class DynamicsKernel<1, 3>;
template class DynamicsKernel<1, 8>;

template class DynamicsKernel<3, 3>;
template class DynamicsKernel<3, 8>;

template class DynamicsKernel<6, 3>;
template class DynamicsKernel<6, 8>;

}
