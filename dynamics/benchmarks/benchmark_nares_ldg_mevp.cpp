/*!
 * @file benchmark_nares_ldg_mevp.cpp
 *
 * Nares channel with the mEVP solver and the LDG discretization of the momentum equation
 * (LDGMEVPDynamicsKernel), as run/config_nares.cfg with
 * DynamicsModule = Nextsim::LDGMEVPDynamics. The velocity space is set by the CMake option
 * LDGVelocity.
 */

#include "LDGMEVPDynamicsKernel.hpp"
#include "benchmark_nares.hpp"

int main(int argc, char* argv[])
{
    using namespace Nextsim;
    const Benchmark::Settings settings
        = Benchmark::parseArguments(argc, argv, "output_nares_ldg_mevp"); // default output directory
    Nares::Problem problem(settings);

    VPParameters params;
    params.nSteps = settings.nSteps;
    params.fc = 0.0; // no Coriolis force
    LDGParameters ldgParams;
    LDGMEVPDynamicsKernel<DGCOMP, LDGCOMP> kernel(params, ldgParams);

    std::cout << "LDG velocity space DG" << LDGCOMP << std::endl;
    Nares::run(kernel, problem, "LDG mEVP");
    return 0;
}
