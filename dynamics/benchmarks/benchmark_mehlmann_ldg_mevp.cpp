/*!
 * @file benchmark_mehlmann_ldg_mevp.cpp
 * @author Thomas Richter <thomas.richter@ovgu.de>
 *
 * Mehlmann benchmark with the mEVP solver and the LDG discretization of the momentum equation
 * (LDGMEVPDynamicsKernel). The velocity space is set by the CMake option LDGVelocity.
 */

#include "LDGMEVPDynamicsKernel.hpp"
#include "benchmark_mehlmann.hpp"

int main(int argc, char* argv[])
{
    using namespace Nextsim;
    const Benchmark::Settings settings
        = Benchmark::parseArguments(argc, argv, "benchmark_mehlmann_ldg_mevp");
    Benchmark::Problem problem(settings);

    VPParameters params;
    params.nSteps = settings.nSteps;
    LDGParameters ldgParams;
    LDGMEVPDynamicsKernel<DGCOMP, LDGCOMP> kernel(params, ldgParams);

    std::cout << "LDG velocity space DG" << LDGCOMP << std::endl;
    Benchmark::run(kernel, problem, "LDG mEVP");
    return 0;
}
