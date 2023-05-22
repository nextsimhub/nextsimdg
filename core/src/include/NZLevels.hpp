/*!
 * @file NZLevels.hpp
 *
 * @date 26 Jan 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef NZLEVELS_HPP
#define NZLEVELS_HPP

#include <cstddef>

namespace Nextsim {

class NZLevels {
public:
    // default constructors and destructors
    static void set(size_t n);
    static size_t get();

private:
    static size_t nZLevels;
};

} /* namespace Nextsim */

#endif /* NZLEVELS_HPP */
