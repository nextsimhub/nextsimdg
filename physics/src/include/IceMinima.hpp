/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ICEMINIMA_HPP
#define ICEMINIMA_HPP

namespace Nextsim {

class IColumnPhysics;
//! A class to hold the minimum ice thresholds without having to pull in
//! IceGrowth and its dependencies.
class IceMinima {
public:
    static inline double h() { return hMin; };
    static inline double c() { return cMin; };
    static const double hMinDefault;
    static const double cMinDefault;

private:
    static double hMin;
    static double cMin;

    friend IColumnPhysics;
};

} /* namespace Nextsim */

#endif /* ICEMINIMA_HPP */
