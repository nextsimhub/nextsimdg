/*----------------------------   testtools.hpp     ---------------------------*/
#ifndef __testtools_HPP
#define __testtools_HPP
/*----------------------------   testtools.hpp     ---------------------------*/

/*!
* @file   testtools.hpp
* @author Thomas Richter <thomas.richter@ovgu.de>
* @date   27/10/2021
*/

#include <array>
#include <vector>

namespace Nextsim {

/*!
   * 
   * Computes the  experimental order of convergence
   *
   * @param values is list of numerical values on meshes (h, h/2, h/4, ...)
   * @return {C,a,q}  a(h) = a + C * h^q
   */
std::array<double, 3> extrapolate(const std::vector<double>& values)
{
    if (values.size() < 3) {
        std::cerr << "TESTTOOLS.HPP 'values' must contain at least 3 numbers" << std::endl;
        abort();
    }
    double a3 = values[values.size() - 1];
    double a2 = values[values.size() - 2];
    double a1 = values[values.size() - 3];

    double D = a1 + a3 - 2.0 * a2;

    if ((D == 0) || (a1 == a2))
        return { 0, 0, -1 };

    double L = (a2 - a3) / (a1 - a2);

    if (L <= 0)
        return { 0, 0, -1 };

    return { (a1 * a1 - 2. * a1 * a2 + a2 * a2) / D, (a1 * a3 - a2 * a2) / D, -log(L) / log(2.0) };
}

}

/*----------------------------   testtools.hpp     ---------------------------*/
/* end of #ifndef __testtools_HPP */
#endif
/*----------------------------   testtools.hpp     ---------------------------*/
