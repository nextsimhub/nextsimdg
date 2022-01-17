/*!
 * @file DevStep.hpp
 *
 * @date Jan 12, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_DEVSTEP_HPP
#define CORE_SRC_INCLUDE_DEVSTEP_HPP

#include "include/Iterator.hpp"

namespace Nextsim {

class DevStep : public Iterator::Iterant {
public:
    DevStep() = default;
    virtual ~DevStep() = default;

    void init(const Environment& env) override {};
    void start(const Iterator::TimePoint& startTime) override {};
    void iterate(const Iterator::Duration& dt) override;
    void stop(const Iterator::TimePoint& stopTime) override {};
};

} /* namespace Nextsim */

#endif /* CORE_SRC_INCLUDE_DEVSTEP_HPP */
