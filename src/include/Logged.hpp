/*
 * Logged.hpp
 *
 *  Created on: 12 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_LOGGED_HPP
#define SRC_INCLUDE_LOGGED_HPP

#include <string>

namespace Nextsim {

class Logged {
public:
	enum level {
		INFO, DEBUG, NOTICE, WARNING, ERROR, CRITICAL, ALERT, EMERGENCY
	};
	static void log(const std::string& message, const level lvl);
	static void info(const std::string& message);
	// TODO: functions for all the levels in between
	static void emergency(const std::string& message);
protected:
	Logged();
	// TODO: Add implementation to actually do some logging
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_LOGGED_HPP */
