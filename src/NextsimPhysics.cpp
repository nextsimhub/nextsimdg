/*!
 * @file NextsimPhysics.cpp
 * @date Sep 9, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/NextsimPhysics.hpp"
#include "include/PrognosticData.hpp"
#include "include/PhysicsData.hpp"

namespace Nextsim {

NextsimPhysics::NextsimPhysics()
{
    // TODO Auto-generated constructor stub

}

NextsimPhysics::~NextsimPhysics()
{
    // TODO Auto-generated destructor stub
}

void NextsimPhysics::updateDerivedData(const PrognosticData& prog, PhysicsData& phys)
{

}

void NextsimPhysics::massFluxOpenWater(const PrognosticData& prog, PhysicsData& phys)
{
    double specificHumidityDifference = phys.specificHumidityWater() - phys.specificHumidityAir();
    phys.evaporationRate() = phys.dragOcean_q() * phys.airDensity() * phys.windSpeed()
            * specificHumidityDifference;
}

void NextsimPhysics::momentumFluxOpenWater(const PrognosticData& prog, PhysicsData& phys)
{
    phys.dragPressure() = phys.airDensity() * phys.dragOcean_m();
}

void NextsimPhysics::heatFluxOpenWater(const PrognosticData& prog, PhysicsData& phys)
{

}
} /* namespace Nextsim */
