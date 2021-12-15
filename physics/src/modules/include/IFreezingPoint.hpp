/*!
 * @file IFreezingPoint.hpp
 *
 * @date Nov 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_IFREEZINGPOINT_HPP_
#define SRC_INCLUDE_IFREEZINGPOINT_HPP_

namespace Nextsim {

//! The interface class for calculation of the freezing point of seawater.
class IFreezingPoint {
public:
    virtual ~IFreezingPoint() = default;

    /*!
     * @brief A virtual function that calculates the freezing point of
     * seawater.
     *
     * @detailed Freezing point in ˚C of water with the given salinity at
     * standard pressure.
     *
     * @param sss Sea surface salinity [PSU]
     */
    virtual double operator()(double sss) const = 0;
};
}
#endif /* SRC_INCLUDE_IFREEZINGPOINT_HPP_ */
