/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ITERATOR_HPP
#define SRC_INCLUDE_ITERATOR_HPP

#include <memory>

#include "Logged.hpp"
#include "Time.hpp"

namespace Nextsim {

//! A class that controls how time steps are performed.
class Iterator {
public:
    typedef TimePoint::Clock Clock;
    class Iterant;

    Iterator() = delete;

    //! Construct a new Iterator given a pointer to an Iterant.
    Iterator(Iterant& iterant)
        : iterant(iterant)
    {
    }

    /*!
     * @brief Sets the time parameters as a start time, stop time and timestep
     * length.
     *
     * @param startTime Start time point.
     * @param stopTime Stop time point.
     * @param timestep Timestep length.
     */
    void setStartStopStep(TimePoint startTime, TimePoint stopTime, Duration timestep);

    void run();

private:
    Iterant& iterant;
    TimePoint startTime;
    TimePoint stopTime;
    Duration timestep;

public:
    //! A base class for classes that specify what happens during one timestep.
    class Iterant {
    public:
        // Define the constructors and copy operator as default to be
        // rule of 5 compliant, given the virtual destructor
        Iterant() = default;
        Iterant(const Iterant& copyFrom) = default;
        Iterant& operator=(const Iterant& copyFrom) = default;
        Iterant(Iterant&&) = default;
        Iterant& operator=(Iterant&&) = default;

        virtual ~Iterant() = default;

        //! Initializes the model, based on some environment stored in the implementing class.
        virtual void init() = 0;
        /*!
         * Initializes the iterant based on the start time.
         *
         * @param startTime the time at the initialization of the iterant.
         */
        virtual void start(const TimePoint& startTime) = 0;
        /*!
         * Performs one iteration a specified length
         *
         * @param dt The length of the timestep.
         */
        virtual void iterate(const TimestepTime& dt) = 0;
        /*!
         * Finalizes the iterant based on the stop time.
         *
         * @param stopTime the time at the finalization of the iterant.
         */
        virtual void stop(const TimePoint& stopTime) = 0;
    };
};
} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ITERATOR_HPP */
