/*!
 * @file ThermoWinton.hpp
 *
 * @date Sep 30, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef THERMOWINTON_HPP
#define THERMOWINTON_HPP

#include "include/Configured.hpp"
#include "include/IIceThermodynamics.hpp"
#include "include/ModelArrayRef.hpp"
namespace Nextsim {

//! A class implementing IIceThermodynamics as the Winton thermodynamics model.
class ThermoWinton : public IIceThermodynamics, public Configured<ThermoWinton> {
public:
    static const size_t nLevels;

    ThermoWinton();
    virtual ~ThermoWinton() = default;

    enum {
        KS_KEY,
        I0_KEY
    };
    void configure() override;

    ModelState getStateRecursive(const OutputSpec& os) const override;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void setData(const ModelState::DataMap&) override;
    void update(const TimestepTime& tsTime) override;

private:
    void calculateElement(size_t i, const TimestepTime& tst);

    HField snowMelt;
    HField topMelt;
    HField botMelt;
    HField qic;
    ModelArrayRef<ProtectedArray::HTRUE_ICE, MARConstBackingStore> oldHi;
    ModelArrayRef<ProtectedArray::SW_IN, MARConstBackingStore> sw_in;
    ModelArrayRef<SharedArray::SUBLIM, MARBackingStore, RO> subl;

    static double i0;
    static const double cVol;

    void calculateTemps(double& tSurf, double& tMidt, double& tBotn, double& mSurf, size_t i, double dt);

};

} /* namespace Nextsim */

#endif /* THERMOWINTON_HPP */
