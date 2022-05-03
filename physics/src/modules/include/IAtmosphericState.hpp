/*!
 * @file IAtmosphericState.hpp
 *
 * @date May 2, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IATMOSPHERICSTATE_HPP
#define IATMOSPHERICSTATE_HPP

#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class IAtmosphericState : public ModelComponent {
public:
    virtual ~IAtmosphericState() = default;

    std::string getName() const override { return "IAtmosphericState"; }
    void setData(const ModelState&) override { }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override {return getState(); }
    std::set<std::string> hFields() override{ return {"sh_air", "sh_water", "sh_ice", "rho_air", "p_air", "cp_air"}; }

    virtual void update(const TimestepTime&) = 0;
protected:
    IAtmosphericState()
    {
        registerModule();
    }

    // Owned, shared fields
    HField sh_air;
    HField sh_water;
    HField sh_ice;
    HField rho_air;
    HField cp_wet;

    ModelArrayRef<ProtectedArray::T_AIR> t_air;
    ModelArrayRef<ProtectedArray::DEW_2M> t_dew2;
    ModelArrayRef<ProtectedArray::SST> sst;
    ModelArrayRef<ProtectedArray::SSS> sss;
    ModelArrayRef<ProtectedArray::T_ICE> tice;
    ModelArrayRef<ProtectedArray::P_AIR> p_air;
};

}

#endif /* IATMOSPHERICSTATE_HPP */
