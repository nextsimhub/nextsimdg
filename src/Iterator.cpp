/*
 * Iterator.cpp
 *
 *  Created on: 11 Aug 2021
 *      Author: Tim Spain, <timothy.spain@nersc.no>
 */

#include "Iterator.hpp"

namespace Nextsim {

Iterator::NullIterant Iterator::nullIterant;

Iterator::Iterator() :
	iterant(&nullIterant)
{ }

Iterator::Iterator(Iterant* iterant) :
		iterant(iterant)
{ }

void Iterator::setIterant(Iterant* iterant) {
	this->iterant = iterant;
}

void Iterator::setStartStopStep(Iterator::TimePoint startTime,
		Iterator::TimePoint stopTime,
		Iterator::Duration timestep) {
	this->startTime = startTime;
	this->stopTime = stopTime;
	this->timestep = timestep;
}

void Iterator::run() {
	iterant->start(startTime);

	for (auto t = startTime; t < stopTime; t += timestep) {
		iterant->iterate(timestep);
	}

	iterant->stop(stopTime);
}

} /* namespace Nextsim */
