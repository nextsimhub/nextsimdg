/*!
 * @file FakeSmeshData.hpp
 *
 * @date Dec 14, 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef FAKESMESHDATA_HPP
#define FAKESMESHDATA_HPP

#include "include/ModelState.hpp"

#include <vector>

namespace Nextsim {

//! A class providing ModelArrays equivalent to those obtained from a restart file.
class FakeSmeshData {
public:
    static ModelState getData();

private:
    static std::vector<double> landmask25km_NH();
};

}

#endif /* FAKESMESHDATA_HPP */
