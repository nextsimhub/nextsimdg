//
// Created by Einar Ólason on 01/09/2022.
//

#ifndef MU71ATMOSPHERE_HPP
#define MU71ATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"
#include "include/IIceAlbedo.hpp"
#include "include/MonthlyCubicBSpline.hpp"
#include "include/constants.hpp"

namespace Nextsim {

/*!
 * @brief The implementation class for the a seasonal atmospheric fluxes following Maykut and
 * Untersteiener's (1971) table 1. Only useful for comparison with that paper and derived setups.
 */
class MU71Atmosphere : public IAtmosphereBoundary, public Configured<MU71Atmosphere> {

public:
    MU71Atmosphere();

    /*!
     * @brief The required update call for an IFluxCalculation implementation. Here we just call
     * calculateElement inside an overElements loop.
     * @param tst The TimestepTime object for the current time step
     */
    void update(const TimestepTime& tst) override;

    enum {
        I0_KEY,
    };
    void configure() override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

private:
    /*!
     * @brief A function to calculate the fluxes qia, qio, and qow (=0), and subl (=0), as well as
     * dqia_dt. All incoming fluxes are tabulated from Maykut and Untersteiner (1971). Outgoing long
     * wave and the derivative (dqio_dt) are black body radiation.
     * @param i The index of the current grid cell
     * @param tst The TimestepTime object for the current time step
     */
    void calculateElement(size_t i, const TimestepTime& tst);

    ModelArrayRef<Protected::T_ICE> tice;
    ModelArrayRef<Protected::HTRUE_SNOW> h_snow_true; // cell-averaged value

    /*!
     * @brief A function to calculate the snow fall according tu Maykut and Untersteiner (1971)
     * @return Snow fall in m/s
     */
    double snowfall();

    // Monthly fluxes from Maykut and Untersteiner (1971)
    const std::vector<double> swTable
        = { 0.00, 0.00, 1.90, 9.99, 17.7, 19.2, 13.6, 9.00, 3.70, 0.40, 0.00, 0.00 };
    const std::vector<double> lwTable
        = { 10.4, 10.3, 10.3, 11.6, 15.1, 18.0, 19.1, 18.7, 16.5, 13.9, 11.2, 10.9 };
    const std::vector<double> shTable
        = { 1.18, 0.76, 0.72, 0.29, -.45, -.39, -.30, -.40, -.17, 0.10, 0.56, 0.79 };
    const std::vector<double> lhTable
        = { 0.00, -.02, -.03, -.09, -.46, -.70, -.64, -.66, -.39, -.19, -.01, -.01 };
    //       Jan,  Feb,  Mar,  Apr,  Mai,  Jun,  Jul,  Aug,  Sept, Oct,  Nov,  Dec

    // Conversion factor from kcal/cm^2/month to W/m^2
    const double convFactor = 4.184e7 / (365.2425 / 12. * 24. * 3600.);

    double dayOfYear;
    bool isLeap;

    monthlyCubicBSpline q_sw;
    monthlyCubicBSpline q_lw;
    monthlyCubicBSpline q_sh;
    monthlyCubicBSpline q_lh;

    static double m_I0;

    IIceAlbedo* iIceAlbedoImpl;
};

}

#endif // NEXTSIM_DG_MU71ATMOSPHERE_HPP
