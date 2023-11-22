/*!
 * @file IceMinima.hpp
 *
 * @date 8 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICEMINIMA_HPP
#define ICEMINIMA_HPP

namespace Nextsim {

class IceGrowth;
//! A class to hold the minimum ice thresholds without having to pull in
//! IceGrowth and its dependencies.
class IceMinima {
public:
    static inline double h() { return hMin; };
    static inline double c() { return cMin; };

private:
    static double hMin;
    static double cMin;
    static const double hMinDefault;
    static const double cMinDefault;

    friend IceGrowth;
};

} /* namespace Nextsim */

#endif /* ICEMINIMA_HPP */
