/*!
 * @file ThermoWinton.cpp
 *
 * @date Sep 30, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoWinton.hpp"

namespace Nextsim {

const size_t ThermoWinton::nLevels = 3;

ThermoWinton::ThermoWinton() { }

void ThermoWinton::setData(const ModelState::DataMap& state)
{
    // The Winton scheme requires three temperature levels in the ice
    if (tice0.data().size() != nLevels * hice.data().size()) {
        double actualLevels = static_cast<double>(tice0.data().size()) / hice.data().size();
        throw std::length_error(std::string("The inferred number of ice temperature levels is ")
            + std::to_string(actualLevels) + " when the Winton ice thermodynamics scheme expects "
            + std::to_string(nLevels));
    }
}

void ThermoWinton::update(const TimestepTime& tst)
{
    overElements(std::bind(&ThermoWinton::calculateElement, this, std::placeholders::_1,
                     std::placeholders::_2),
        tst);
}

void ThermoWinton::calculateElement(size_t i, const TimestepTime& tst)
{
    // Calculate conductivities and other coefficients
    CoeffArray c;
    calculateCoeffs(c, i);

    // Is the surface melting?

    // Recalculate T1 and Tsurf if so

    // Calculate T2

    // Thickness changes
    // ice

    // snow

    // sublimation
    // 4 cases

    // Bottom melt/freezing

    // Melting at the surface

    // Snow to ice conversion

    // Adjust the temperatures to evenly divide the ice

    // Remove very small ice thickness
}

void ThermoWinton::calculateCoeffs(CoeffArray& c, size_t i)
{
    c[K12] = 0; // &c.
    c[A] = qia[i] - tice0.zIndexAndLayer(i, 0) * qQia_dt[i];
    c[B] = dQia_dt[i];
}

} /* namespace Nextsim */
