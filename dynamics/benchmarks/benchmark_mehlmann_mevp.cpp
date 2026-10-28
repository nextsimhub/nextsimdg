/*!
 * @file benchmark_mehlmann_mevp.cpp
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 * Mehlmann benchmark with the mEVP solver and CG velocities (MEVPDynamicsKernel).
 */

#include "MEVPDynamicsKernel.hpp"
#include "benchmark_mehlmann.hpp"

int main(int argc, char* argv[])
{
    using namespace Nextsim;
    const Benchmark::Settings settings
        = Benchmark::parseArguments(argc, argv, "benchmark_mehlmann_mevp");
    Benchmark::Problem problem(settings);

    VPParameters params;
    params.nSteps = settings.nSteps;
    MEVPDynamicsKernel<DGCOMP> kernel(params);

    Benchmark::run(kernel, problem, "CG mEVP");
    return 0;
}
