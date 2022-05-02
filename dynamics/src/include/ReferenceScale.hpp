/*!
 * @file ReferenceScale.hpp
 * @date 1 Mar 2022
 * @author Thomas Richter <thomas.richter@ovgu.de>
 */

#ifndef __REFERENCESCALE_HPP
#define __REFERENCESCALE_HPP

namespace Nextsim {

//! Reference scaling
class ReferenceScale {
public:
    //! Spatial scale. x=1 [no unit] refers to L [m]
    static constexpr double L = 512.e3; // 1  |-> 512 km

    //! Temporal scale. t=1 [no unit] refers to T [sec]
    static constexpr double T = 1.e3; // 1  |-> 1000 sec approx 16 min.

    //! Stress scale. S=1 [no unit] refers to S [Pa] = [kg/m/s^2]
    static constexpr double S = 1.e5; // 1  |-> 100 kPa
};

} /* namespace Nextsim */

#endif /* __REFERENCESCALE_HPP */
