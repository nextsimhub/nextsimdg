/*!
 * @file Time.hpp
 *
 * @date Mar 15, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_TIME_HPP
#define CORE_SRC_INCLUDE_TIME_HPP

#include <chrono>

namespace Nextsim {

typedef int TimePoint; // TODO Use a real time type
typedef int Duration; // TODO Use a real duration type
                      //    typedef std::chrono::time_point<Clock> TimePoint;
                      //    typedef std::chrono::seconds Duration;

}; // namespace Nextsim

#endif /* CORE_SRC_INCLUDE_TIME_HPP */
