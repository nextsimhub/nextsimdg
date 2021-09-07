/**
 * @file Timed.cpp
 * @date 11 Aug 2021
 * @author Tim Spain, <timothy.spain@nersc.no>
 */

#include <string>

#include "Timed.hpp"
#include "Timer.hpp"

namespace Nextsim {
Timer staticTimer = Timer();
Timer& Timed::timer(staticTimer);

Timed::Timed() { }

void Timed::tick(const std::string& timerName) {
	timer.tick(timerName);
}

void Timed::tock(const std::string& timerName) {
	timer.tock(timerName);
}

std::string Timed::report() {
	return timer.report();
}
}
