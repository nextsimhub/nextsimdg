/*
 * Timed_test.cpp
 *
 *  Created on: 11 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN

#include <catch2/catch.hpp>

#include <string>
#include <sstream>

#include "../src/include/Timed.hpp"

namespace Nextsim {

class TimeMe: public Timed {
public:
	TimeMe() { };
	static const std::string timerName;
	void tickTockItsTimingOClock( ) {
		tick(timerName);
		tock(timerName);
	}
};

const std::string TimeMe::timerName = "tickTock";

// Requires use of the Timer defined in CountTimer.cpp
TEST_CASE( "Count timer testing", "[Timed]") {
	TimeMe timeMe;

	std::stringstream builder;
	builder << "Timing:" << std::endl;
	std::string reported = timeMe.report();

	REQUIRE(builder.str() == reported);

	timeMe.tickTockItsTimingOClock();
	reported = timeMe.report();

	// Target string taken from CountTimer.cpp:report()
	auto key = TimeMe::timerName;
	builder << key << " started " << 1 << ", stopped " << 1 << " times" << std::endl;

	REQUIRE(builder.str() == reported);
}

} /* namespace Nextsim */
