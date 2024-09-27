/*!
 * @file dgLimit.hpp
 * @date 27 September 2024
 * @author Robert Jendersie <robert.jendersie@ovgu.de>
 */

#ifndef __KOKKOSDGLIMIT_HPP
#define __KOKKOSDGLIMIT_HPP

#include "KokkosUtils.hpp"
#include "../../include/dgVector.hpp"

namespace Nextsim {

void limitMax(const KokkosDeviceView<DGVector<1>>& dg, FloatType max);
void limitMax(const KokkosDeviceView<DGVector<3>>& dg, FloatType max);
void limitMax(const KokkosDeviceView<DGVector<6>>& dg, FloatType max);

void limitMin(const KokkosDeviceView<DGVector<1>>& dg, FloatType min);
void limitMin(const KokkosDeviceView<DGVector<3>>& dg, FloatType min);
void limitMin(const KokkosDeviceView<DGVector<6>>& dg, FloatType min);

}

#endif // __KOKKOSDGLIMIT_HPP