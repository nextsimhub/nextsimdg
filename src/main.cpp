/**
 * @file main.cpp
 * @date 11 Aug 2021
 * @author Tim Spain, <timothy.spain@nersc.no>
 */

#include <iostream>

#include "Model.hpp"

int main(int argc, char* argv[]) {

	// Construct the Model
	Nextsim::Model model = Nextsim::Model();
	// Run the Model
	model.run();

	return 0;
}


