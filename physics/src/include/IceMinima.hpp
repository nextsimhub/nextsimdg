/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICEMINIMA_HPP
#define ICEMINIMA_HPP

#include "include/FloatType.hpp"

namespace Nextsim {

class IColumnPhysics;
//! A class to hold the minimum ice thresholds without having to pull in
//! IColumnPhysics and its dependencies.
class IceMinima {
public:
    static inline FloatType h() { return hMin; };
    static inline FloatType c() { return cMin; };
    static const FloatType hMinDefault;
    static const FloatType cMinDefault;

private:
    static FloatType hMin;
    static FloatType cMin;

    friend IColumnPhysics;
};

} /* namespace Nextsim */

#endif /* ICEMINIMA_HPP */
