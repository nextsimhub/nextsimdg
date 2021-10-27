/*!
 * @file TimerBase.hpp
 * @date 25 Oct 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_TIMERBASE_HPP
#define SRC_INCLUDE_TIMERBASE_HPP

#include <ostream>
#include <string>

namespace Nextsim {

class TimerBase {
public:
    typedef std::string Key;

    virtual ~TimerBase() = default;

    virtual void tick(const Key& timerName) = 0;
    virtual void tock(const Key& timerName) = 0;

    virtual double lap(const Key& timerName) const = 0;
    virtual double elapsed(const Key& timerName) const = 0;

    virtual std::ostream& report(const Key& timerName, std::ostream& os) const = 0;
    virtual std::ostream& report(std::ostream& os) const = 0;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_TIMERBASE_HPP */
