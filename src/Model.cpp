 /**
  * @file Model.cpp
  * @date 12 Aug 2021
  * @author Tim Spain, <timothy.spain@nersc.no>
  */
	
#include "Model.hpp"

#include "SimpleIterant.hpp"

namespace Nextsim {

Model::Model( )
{
	iterant = new SimpleIterant( );
	deleteIterant = true;
	iterator.setIterant(iterant);

	const int runLength = 5;

	Iterator::TimePoint now(std::chrono::system_clock::now());
	Iterator::Duration dt = std::chrono::seconds(1);
	Iterator::TimePoint hence = now + runLength * dt;

	iterator.setStartStopStep(now, hence, dt);
}

// TODO: add another constructor which takes arguments specifying the
// environment and configuration. This will be the location of the
// logic with selects the components of the model that will run,
// translates I/O details from file configuration to object variable
// values, specifies file paths and likely more besides.

Model::~Model( ) {
	if (deleteIterant) {
		delete iterant;
	}
}

void Model::run( ) {
	iterator.run();
}
} /* namespace Nextsim */
