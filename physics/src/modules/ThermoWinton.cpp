/*!
 * @file ThermoWinton.cpp
 *
 * @date Sep 30, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ThermoWinton.hpp"

#include "include/constants.hpp"

#include <cmath>

namespace Nextsim {

const size_t ThermoWinton::nLevels = 3;
double ThermoWinton::IIceThermodynamics::kappa_s;
double ThermoWinton::i0;
static const double k_sDefault = 0.3096;
static const double i0_default = 0.17;

const double ThermoWinton::cVol = Ice::cp * Ice::rho; // bulk heat capacity of ice

ThermoWinton::ThermoWinton() { }

template <>
const std::map<int, std::string> Configured<ThermoWinton>::keyMap = {
    { ThermoWinton::KS_KEY, IIceThermodynamics::getKappaSConfigKey() },
    { ThermoWinton::I0_KEY, "nextsim_thermo.I_0" },
};

void ThermoWinton::configure()
{
    kappa_s = Configured::getConfiguration(keyMap.at(KS_KEY), k_sDefault);
    i0 = Configured::getConfiguration(keyMap.at(I0_KEY), i0_default);
}

ModelState ThermoWinton::getStateRecursive(const OutputSpec& os) const
{
    ModelState state = { {},
        {
            { keyMap.at(KS_KEY), kappa_s },
            { keyMap.at(I0_KEY), i0 },
        } };
    return os ? state : ModelState();
}

ThermoWinton::HelpMap& ThermoWinton::getHelpText(HelpMap& map, bool getAll)
{
    map["ThermoWinton"] = {
        { keyMap.at(KS_KEY), ConfigType::NUMERIC, { "0", "∞" }, std::to_string(k_sDefault),
            "W K⁻¹ m⁻¹", "Thermal conductivity of snow." },
        { keyMap.at(I0_KEY), ConfigType::NUMERIC, { "0", "1" }, std::to_string(i0_default),
            "unitless", "Optical albedo of liquid water." },
        };
    return map;
}
ThermoWinton::HelpMap& ThermoWinton::getHelpRecursive(HelpMap& map, bool getAll)
{
    return getHelpText(map, getAll);
}

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
    double tSurf = tice0.zIndexAndLayer(i, 0);
    double tMidt = tice0.zIndexAndLayer(i, 1);
    double tBotn = tice0.zIndexAndLayer(i, 2);

    double mSurf = 0; // surface melting mass loss
    // Calculate temperatures by solving the heat conduction equation
    calculateTemps(tSurf, tMidt, tBotn, mSurf, i, tst.step.seconds());

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

void ThermoWinton::calculateTemps(double& tSurf, double& tMidt, double& tBotn, double& mSurf, size_t i, double dt)
{

    double& hi = hice[i];

    double tMelt = (hsnow[i] > 0) ? 0 : Ice::Tm; // Melting point at the surface

    double k12
        = 4 * Ice::kappa * kappa_s / (kappa_s * hi + 4 * Ice::kappa * hsnow[i]); // Winton & al. (5)
    double a = qia[i] - tice0.zIndexAndLayer(i, 0) * dQia_dt[i];
    double b = dQia_dt[i];
    double k32 = 2 * Ice::kappa / hi;

    double a1 =
            hi * cVol / (2 * dt) +
            k32 * (4 * dt * k32 + hi * cVol) / (6 * dt * k32 + hi * cVol) +
            b * k12 / (k12 + b);
    double b1 =
            -hi * (cVol * tMidt + Ice::Lf * Ice::rho / tMidt) / (2 * dt) -
            i0 * sw_in[i] -
            k32 * (4 * dt * k32 * tf[i] + hi * cVol * tBotn) / (6 * dt * k32 + hi * cVol) +
            a * k12 / (k12 +b);
    double c1 = hi * Ice::Lf * Ice::rho * Ice::Tm / (2 * dt);

    // Updated surface and mid-ice temperatures
    tMidt = -(b1 + std::sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1);
    tSurf = (k12 * tMidt - a) / (k12 + b);

    // Is the surface melting?
    if (tSurf > tMelt) {
        // Recalculate T1 and Tsurf if so
        tSurf = tMelt;
        // apply the change to the *1 parameters, rather than recalculating in full
        a1 += k12 - k12 * b / (k12 + b);
        b1 -= k12 * tSurf + a * k12 / (k12 + b);
        tMidt = -(b1 + std::sqrt(b1 * b1 - 4 * a1 * c1)) / (2 * a1);

        mSurf = k12 * (tMidt - tSurf) - (a + b * tSurf);
    }

    // update T_botn based on the new value of tMidt
    tBotn = (2 * dt * k32 * (tMidt + 2 * tf[i]) + hi * cVol * tBotn) / (6 * dt * k32 + hi * cVol);
}

} /* namespace Nextsim */
