/*!
 * @file Iterator.hpp
 * @date 11 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_ITERATOR_HPP
#define SRC_INCLUDE_ITERATOR_HPP

#include "Logged.hpp"
#include "Time.hpp"

namespace Nextsim {

//! A class that controls how time steps are performed.
class Iterator {
public:
    typedef TimePoint::Clock Clock;
    class Iterant;

    Iterator();
    //! Construct a new Iterator given a pointer to an Iterant.
    Iterator(Iterant* iterant);

    /*!
     * @brief Sets the iterant to be iterated using a pointer.
     *
     * @param iterant The Iterant that defines a single timestep of the model.
     */
    void setIterant(Iterant* iterant);

    /*!
     * @brief Sets the time parameters as a start time, stop time and timestep
     * length.
     *
     * @param startTime Start time point.
     * @param stopTime Stop time point.
     * @param timestep Timestep length.
     */
    void setStartStopStep(TimePoint startTime, TimePoint stopTime, Duration timestep);
    /*!
     * @brief Sets the time parameters as a start time, run length and timestep
     * length.
     *
     * @param startTime Start time point.
     * @param duration Minimum length of the run.
     * @param timestep Timestep length.
     */
    void setStartDurationStep(TimePoint startTime, Duration duration, Duration timestep);

    /*!
     * @brief Parses the four strings and sets the time parameters from them.
     *
     * @details Give the four strings, the parser determines which ones are
     * valid and can be used. The strings are then parsed to time values and
     * the start time returned, as this can be considered the current time.
     *
     * @param startTimeStr string to parse for the model start time.
     * @param stopTimeStr string to parse for the model stop time.
     * @param durationStr string to parse for the model run duration.
     * @param stepStr string to parse for the model time step length.
     */
    TimePoint parseAndSet(const std::string& startTimeStr, const std::string& stopTimeStr,
        const std::string& durationStr, const std::string& stepStr);
    //! Run the Iterant over the specified time period.
    void run();

private:
    Iterant* iterant;
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

    //! A simple Iterant that does nothing.
    class NullIterant : public Iterant {
        inline void init() {};
        inline void start(const TimePoint& startTime) {};
        inline void iterate(const TimestepTime& dt) {};
        inline void stop(const TimePoint& stopTime) {};
    };
};
} /* namespace Nextsim */

#endif /* SRC_INCLUDE_ITERATOR_HPP */
