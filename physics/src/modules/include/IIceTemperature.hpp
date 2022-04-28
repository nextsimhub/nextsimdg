/*!
 * @file IIceTemperature.hpp
 *
 * @date Apr 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IICETEMPERATURE_HPP
#define IICETEMPERATURE_HPP

#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class IIceTemperature : public ModelComponent {
public:
    IIceTemperature()
    {
        registerSharedArray(SharedArray::T_ICE, &tice);
        registerSharedArray(SharedArray::Q_IC, &qic);
    }
    virtual ~IIceTemperature() = default;

    void setData(const ModelState&) override {};
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

    std::string getName() const override { return "IIceTemperature"; }

    std::set<std::string> hFields() const override { return { "Qic" }; }
    std::set<std::string> zFields() const override { return { "Tice" }; }

    virtual void update(const TimestepTime& tst) = 0;

protected:
    // Owned shared arrays
    ZField tice;
    HField qic;

    // Referenced arrays
    // Prognostic
    ModelArrayRef<ProtectedArray::T_ICE> tice0;
    ModelArrayRef<ProtectedArray::TF> tf;
    ModelArrayRef<ProtectedArray::H_ICE> hice;
    ModelArrayRef<ProtectedArray::H_SNOW> hsnow;
    // Derived data
    ModelArrayRef<SharedArray::Q_IA, RO> qia;
    ModelArrayRef<SharedArray::DQIA_DT, RO> dqia_dt;
};

class ConstantIceTemperature : public IIceTemperature {
    void update(const TimestepTime& tst) override { }
};

}

#endif /* IICETEMPERATURE_HPP */
