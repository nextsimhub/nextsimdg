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
    void init();
    //! Prints an informative message about starting.
    void start(const TimePoint& startTime);
    //! Prints an informative message about iterating.
    void iterate(const Duration& dt);
    //! Prints an informative message about stopping.
    void stop(const TimePoint& stopTime);

    template <typename T> static T zeroTime();

private:
    // Functions to deal with different realizations of the Iterator time types
    static int count(const std::chrono::seconds& dt) { return dt.count(); };
    static int count(const int dt) { return dt; };

    static std::string stringFromTimePoint(
        const std::chrono::time_point<std::chrono::system_clock>& t);
    static std::string stringFromTimePoint(const int t);
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_SIMPLEITERANT_HPP */
