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
class ThermoWinton : public IIceThermodynamics {
public:
    static const size_t nLevels;

    ThermoWinton();
    virtual ~ThermoWinton() = default;

    enum {
        KS_KEY,
    };
    void configure() override;

    ModelState getStateRecursive(const OutputSpec& os) const override;

    void setData(const ModelState::DataMap&) override;
    void update(const TimestepTime& tsTime) override;

private:
    void calculateElement(size_t i, const TimestepTime& tst);

    HField snowMelt;
    HField topMelt;
    HField botMelt;
    HField qic;
    ModelArrayRef<ProtectedArray::HTRUE_ICE, MARConstBackingStore> oldHi;

    enum {
        K12, A, B, K32, A1, B1, C1, T1, TSURF, COUNT
    };
    typedef std::array<double, COUNT> CoeffArray;

    static double kappa_s;

    void calculateCoeffs(CoeffArray& c, size_t i, double dt);

};

} /* namespace Nextsim */

#endif /* THERMOWINTON_HPP */
