/*!
 * @file benchmark_mehlmann_bbm.cpp
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 * Mehlmann benchmark with the brittle Bingham-Maxwell rheology and CG velocities
 * (BBMDynamicsKernel).
 */

#include "BBMDynamicsKernel.hpp"
#include "benchmark_mehlmann.hpp"

int main(int argc, char* argv[])
{
    using namespace Nextsim;
    const Benchmark::Settings settings
        = Benchmark::parseArguments(argc, argv, "benchmark_mehlmann_bbm");
    Benchmark::Problem problem(settings);

    BBMParameters params;
    params.nSteps = settings.nSteps;
    BBMDynamicsKernel<DGCOMP> kernel(params);

    Benchmark::run(kernel, problem, "CG BBM", true);
    return 0;
}
