/*!
 * @file Iterator.hpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ITERATOR_HPP
#define SRC_INCLUDE_ITERATOR_HPP

#include <chrono>

#include "Logged.hpp"
#include "Timed.hpp"

namespace Nextsim {

class Environment;

class Iterator : public Timed, public Logged {
public:
    typedef std::chrono::system_clock Clock;
    typedef std::chrono::time_point<Clock> TimePoint;
    typedef std::chrono::seconds Duration;
    class Iterant;

    Iterator();
    Iterator(Iterant* iterant);

    void setIterant(Iterant* iterant);

    void setStartStopStep(TimePoint startTime, TimePoint stopTime, Duration timestep);
    void setStartDurationStep(TimePoint startTime, Duration duration, Duration timestep);

    void run();

private:
    Iterant* iterant; // FIXME smart pointer
    TimePoint startTime;
    TimePoint stopTime;
    Duration timestep;

public:
    class Iterant : public Logged, public Timed {
    public:
        // Define the constructors and copy operator as default to be
        // rule of 5 compliant, given the virtual destructor
        Iterant() = default;
        Iterant(const Iterant& copyFrom) = default;
        Iterant& operator=(const Iterant& copyFrom) = default;
        Iterant(Iterant&&) = default;
        Iterant& operator=(Iterant&&) = default;

        virtual ~Iterant() = default;

        virtual void init(const Environment&) = 0;
        virtual void start(const TimePoint& startTime) = 0;
        virtual void iterate(const Duration& dt) = 0;
        virtual void stop(const TimePoint& stopTime) = 0;
    };

    class NullIterant : public Iterant {
        inline void init(const Environment& env) {};
        inline void start(const Iterator::TimePoint& startTime) {};
        inline void iterate(const Iterator::Duration& dt) {};
        inline void stop(const Iterator::TimePoint& stopTime) {};
    };

    static NullIterant nullIterant;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ITERATOR_HPP */
