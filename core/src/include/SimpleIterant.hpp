/*!
 * @file SimpleIterant.hpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_SIMPLEITERANT_HPP
#define SRC_INCLUDE_SIMPLEITERANT_HPP

#include "Iterator.hpp"

namespace Nextsim {

//! A simple implementation of the Iterator::Iterant class.
class SimpleIterant : public Iterator::Iterant {
public:
    SimpleIterant();

    //! Prints an informative message about initialization.
    void init(const Environment&);
    //! Prints an informative message about starting.
    void start(const Iterator::TimePoint& startTime);
    //! Prints an informative message about iterating.
    void iterate(const Iterator::Duration& dt);
    //! Prints an informative message about stopping.
    void stop(const Iterator::TimePoint& stopTime);
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_SIMPLEITERANT_HPP */
