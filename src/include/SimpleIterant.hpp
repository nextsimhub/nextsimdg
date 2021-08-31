/*
 * SimpleIterant.hpp
 *
 *  Created on: 12 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_SIMPLEITERANT_HPP
#define SRC_INCLUDE_SIMPLEITERANT_HPP

#include "Iterator.hpp"

namespace Nextsim {

class SimpleIterant: public Iterator::Iterant {
public:
	SimpleIterant();

	void init(const Environment&);
	void start(const Iterator::TimePoint& startTime);
	void iterate(const Iterator::Duration & dt);
	void stop(const Iterator::TimePoint& stopTime);

};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_SIMPLEITERANT_HPP */
