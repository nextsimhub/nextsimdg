/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/UniformOcean.hpp"

#include "include/IFreezingPoint.hpp"
#include "include/ModelComponent.hpp"
#include "include/NextsimModule.hpp"
#include "include/constants.hpp"

namespace Nextsim {

TEST_SUITE_BEGIN("UniformOcean");
TEST_CASE("UniformOcean construction")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    FloatType sstIn = -1.23;
    FloatType sssIn = 34.5;
    FloatType mldIn = 12.3;
    FloatType uIn = -0.12;
    FloatType vIn = 0.98;
    FloatType qioIn = 123.45;

    UniformOcean uniOcn(sstIn, sssIn, mldIn, uIn, vIn, qioIn);
    uniOcn.setData(ModelState::DataMap());

    ModelArrayAccessor<Protected::SST> sstAccessor(ModelComponent::getStore());
    const HField& sst = sstAccessor.getHostRO();
    ModelArrayAccessor<Protected::SSS> sssAccessor(ModelComponent::getStore());
    const HField& sss = sssAccessor.getHostRO();
    ModelArrayAccessor<Protected::MLD> mldAccessor(ModelComponent::getStore());
    const HField& mld = mldAccessor.getHostRO();
    ModelArrayAccessor<Protected::OCEAN_U> uAccessor(ModelComponent::getStore());
    const HField& u = uAccessor.getHostRO();
    ModelArrayAccessor<Protected::OCEAN_V> vAccessor(ModelComponent::getStore());
    const HField& v = vAccessor.getHostRO();
    ModelArrayAccessor<Shared::Q_IO, RO> qioAccessor(ModelComponent::getStore());
    const HField& qio = qioAccessor.getHostRO();
    ModelArrayAccessor<Protected::ML_BULK_CP> cpmlAccessor(ModelComponent::getStore());
    const HField& cpml = cpmlAccessor.getHostRO();
    ModelArrayAccessor<Protected::TF> tfAccessor(ModelComponent::getStore());
    const HField& tf = tfAccessor.getHostRO();

    REQUIRE(sst[0] == sstIn);
    REQUIRE(sss[0] == sssIn);
    REQUIRE(mld[0] == mldIn);
    REQUIRE(u[0] == uIn);
    REQUIRE(v[0] == vIn);
    REQUIRE(qio[0] == qioIn);
    REQUIRE(cpml[0] == mldIn * Water::rho * Water::cp);
    REQUIRE(tf[0] == Module::getImplementation<IFreezingPoint>()(sssIn));
}

TEST_CASE("UniformOcean set functions")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    Module::setImplementation<IFreezingPoint>("Nextsim::UnescoFreezing");

    FloatType sstIn = -2.34;
    FloatType sssIn = 32.1;
    FloatType mldIn = 23.4;
    FloatType uIn = 0.12;
    FloatType vIn = -1.23;
    FloatType qioIn = 234.5;

    UniformOcean uniOcn;
    uniOcn.setSST(sstIn).setSSS(sssIn).setMLD(mldIn).setU(uIn).setV(vIn).setQio(qioIn);
    uniOcn.setData(ModelState::DataMap());

    ModelArrayAccessor<Protected::SST> sstAccessor(ModelComponent::getStore());
    const HField& sst = sstAccessor.getHostRO();
    ModelArrayAccessor<Protected::SSS> sssAccessor(ModelComponent::getStore());
    const HField& sss = sssAccessor.getHostRO();
    ModelArrayAccessor<Protected::MLD> mldAccessor(ModelComponent::getStore());
    const HField& mld = mldAccessor.getHostRO();
    ModelArrayAccessor<Protected::OCEAN_U> uAccessor(ModelComponent::getStore());
    const HField& u = uAccessor.getHostRO();
    ModelArrayAccessor<Protected::OCEAN_V> vAccessor(ModelComponent::getStore());
    const HField& v = vAccessor.getHostRO();
    ModelArrayAccessor<Shared::Q_IO, RO> qioAccessor(ModelComponent::getStore());
    const HField& qio = qioAccessor.getHostRO();
    ModelArrayAccessor<Protected::ML_BULK_CP> cpmlAccessor(ModelComponent::getStore());
    const HField& cpml = cpmlAccessor.getHostRO();
    ModelArrayAccessor<Protected::TF> tfAccessor(ModelComponent::getStore());
    const HField& tf = tfAccessor.getHostRO();

    REQUIRE(sst[0] == sstIn);
    REQUIRE(sss[0] == sssIn);
    REQUIRE(mld[0] == mldIn);
    REQUIRE(u[0] == uIn);
    REQUIRE(v[0] == vIn);
    REQUIRE(qio[0] == qioIn);
    REQUIRE(cpml[0] == mldIn * Water::rho * Water::cp);
    REQUIRE(tf[0] == Module::getImplementation<IFreezingPoint>()(sssIn));
}
TEST_SUITE_END();

} /* namespace Nextsim */
