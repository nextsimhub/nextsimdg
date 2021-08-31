/*
 * CountTimer.cpp
 *
 *  Created on: 11 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#include "Timer.hpp"

#include <string>
#include <sstream>
#include <map>
#include <utility>

static std::map<std::string, int> startMap;
static std::map<std::string, int> stopMap;

namespace Nextsim {
Timer::Timer( ) { }

void Timer::tick(const std::string& timerName) {
	if (startMap.count(timerName) == 0) {
		startMap[timerName] = 1;
	} else {
		startMap[timerName]++;
	}
}

void Timer::tock(const std::string& timerName) {
	if (stopMap.count(timerName) == 0) {
		stopMap[timerName] = 1;
	} else {
		stopMap[timerName]++;
	}
}

std::string Timer::report() const {
	std::stringstream builder;
	builder << "Timing:" << std::endl;
	for (auto iter = startMap.begin(); iter != startMap.end(); iter++) {
		auto key = iter->first;
		builder << key << " started " << iter->second << ", stopped " << stopMap[key] << " times" << std::endl;
	}
	return builder.str();
}
}

