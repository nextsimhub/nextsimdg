/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/ModelComponent.hpp"
#include "include/ModelArrayAccessor.hpp"

#include <stdexcept>

namespace Nextsim {

// (Ab)use the exception mechanism to inform doctest that things are working correctly internally.
class HappyExcept : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

class Module1 : public ModelComponent {
public:
    Module1() { }
    std::string getName() const override { return "Module1"; }
    void setData(const ModelState::DataMap& st) override
    {
        throw(HappyExcept(std::string("setData for ") + getName()));
    }
};

class ModuleSupplyAndWait : public ModelComponent {
public:
    ModuleSupplyAndWait()
        : hiceAccessor(getStore(), RO)
        , cice_refAccessor(getStore())
    {
    }
    void setData(const ModelState::DataMap& ms) override { hiceAccessor.getHostRW()[0] = hiceData; }
    std::string getName() const override { return "SupplyAndWait"; }
    ModelState getStatePrognostic() const override
    {
        return { {
                     { "hice", hiceAccessor.getHostRO() },
                 },
            {} };
    }

    const double hiceData = 1.2;
    double data() { return hiceAccessor.getHostRO()[0]; }
    double refData() { return cice_refAccessor.getHostRO()[0]; }

private:
    ModelArrayAccessor<Shared::H_ICE_DG, RW> hiceAccessor;
    ModelArrayAccessor<Shared::C_ICE_DG> cice_refAccessor;
};

class ModuleRequestAndSupply : public ModelComponent {
public:
    ModuleRequestAndSupply()
        : ciceAccessor(getStore(), RO)
        , hice_refAccessor(getStore())
    {
    }
    void setData(const ModelState::DataMap& ms) override { ciceAccessor.getHostRW()[0] = ciceData; }
    std::string getName() const override { return "SupplyAndWait"; }
    ModelState getStatePrognostic() const override
    {
        return { {
                     { "cice", ciceAccessor.getHostRO() },
                 },
            {} };
    }

    const double ciceData = 0.6;
    double data() { return ciceAccessor.getHostRO()[0]; }
    double refData() { return hice_refAccessor.getHostRO()[0]; }

private:
    ModelArrayAccessor<Shared::C_ICE_DG, RW> ciceAccessor;
    ModelArrayAccessor<Shared::H_ICE_DG> hice_refAccessor;
};

TEST_SUITE_BEGIN("ModelComponent");
TEST_CASE("Test array registration")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });
    ModuleSupplyAndWait saw;
    ModuleRequestAndSupply ras;

    REQUIRE(ras.data() == saw.refData());
    REQUIRE(saw.data() == ras.refData());
}

class ModuleSemiShared : public ModelComponent {
public:
    ModuleSemiShared()
        : qicAccessor(getStore(), RW)
        , qio_refAccessor(getStore())
    {
    }
    void setData(const ModelState::DataMap& ms) override { qicAccessor.getHostRW()[0] = qicData; }
    std::string getName() const override { return "SemiShared"; }
    ModelState getStatePrognostic() const override
    {
        return { {
                     { "qic", qicAccessor.getHostRO() },
                 },
            {} };
    }

    const double qicData = 123;
    double data() { return qicAccessor.getHostRO()[0]; }
    double refData() { return qio_refAccessor.getHostRO()[0]; }

private:
    ModelArrayAccessor<Shared::Q_IC, RW> qicAccessor;
    ModelArrayAccessor<Shared::Q_IO, RO> qio_refAccessor;
};

class ModuleShared : public ModelComponent {
public:
    ModuleShared()
        : qioAccessor(getStore(), RW)
        , qic_refAccessor(getStore())
    {
    }
    void setData(const ModelState::DataMap& ms) override { /*qio[0]; */ }
    std::string getName() const override { return "Shared"; }
    ModelState getStatePrognostic() const override
    {
        return { {
                     { "qio", qioAccessor.getHostRO() },
                 },
            {} };
    }

    const double qioData = 234;
    const double qicData = 246;
    double data() { return qioAccessor.getHostRO()[0]; }
    double& refData() { return qic_refAccessor.getHostRW()[0]; }
    void setRefData() { qic_refAccessor.getHostRW()[0] = qicData; }

private:
    ModelArrayAccessor<Shared::Q_IO, RW> qioAccessor;
    ModelArrayAccessor<Shared::Q_IC, RW> qic_refAccessor;
};

TEST_CASE("Shared and semi-protected arrays")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1, 1 });

    ModuleSemiShared semi;
    ModuleShared share;

    REQUIRE(share.data() == semi.refData());
    REQUIRE(semi.data() == share.refData());

    share.refData() = share.qicData;

    REQUIRE(semi.data() == share.qicData);
}
TEST_SUITE_END();

} /* namespace Nextsim */
