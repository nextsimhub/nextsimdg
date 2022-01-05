/*!
 * @file Physics1dBase.hpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_PHYSICS1DBASE_HPP
#define SRC_INCLUDE_PHYSICS1DBASE_HPP

#include "include/ElementData.hpp"

namespace Nextsim {

/*!
 * @brief Base class for 1d column physics.
 *
 * @details The functions in this class define an interface which allows an
 * implementable set of methods to calculate the physics (mass flow, drag,
 * thermodynamics) in a single grid or mesh cell.
 */
class Physics1dBase {
public:
    Physics1dBase() = default;
    virtual ~Physics1dBase() = default;

    //! Performs the one dimensional physics calculations.
    void physics1d(ElementData&);
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_PHYSICS1DBASE_HPP */
