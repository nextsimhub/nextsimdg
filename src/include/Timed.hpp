/*
 * Timed.h
 *
 *  Created on: 11 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_TIMED_HPP
#define SRC_INCLUDE_TIMED_HPP

#include <string>

#include "Timer.hpp"

namespace Nextsim {

class Timed {
public:
	static void tick(const std::string& timerName);
	static void tock(const std::string& timerName);
	static std::string report();
protected:
	Timed();
private:
	static Timer& timer;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_TIMED_HPP */
