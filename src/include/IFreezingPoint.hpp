/*!
 * @file IFreezingPoint.hpp
 *
 * @date Nov 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_IFREEZINGPOINT_HPP_
#define SRC_INCLUDE_IFREEZINGPOINT_HPP_

class IFreezingPoint {
    virtual ~IFreezingPoint() = default;

    //! Freezing point in ˚C of water with salinity sss and standard pressure
    virtual double operator()(double sss) const = 0;
};

#endif /* SRC_INCLUDE_IFREEZINGPOINT_HPP_ */
