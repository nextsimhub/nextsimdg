/*!
 * @file benchmark_nares_mevp.cpp
 *
 * Nares channel with the mEVP solver and CG velocities (MEVPDynamicsKernel), as
 * run/config_nares.cfg with DynamicsModule = Nextsim::MEVPDynamics.
 */

#include "MEVPDynamicsKernel.hpp"
#include "benchmark_nares.hpp"

int main(int argc, char* argv[])
{
    using namespace Nextsim;
    const Benchmark::Settings settings
        = Benchmark::parseArguments(argc, argv, "output_nares_mevp"); // default output directory
    Nares::Problem problem(settings);

    VPParameters params;
    params.nSteps = settings.nSteps;
    params.fc = 0.0; // no Coriolis force
    MEVPDynamicsKernel<DGCOMP> kernel(params);

    Nares::run(kernel, problem, "CG mEVP");
    return 0;
}
