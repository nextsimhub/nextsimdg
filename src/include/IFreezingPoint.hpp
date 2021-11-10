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

    virtual double operator()(double sss) = 0;
};

#endif /* SRC_INCLUDE_IFREEZINGPOINT_HPP_ */
