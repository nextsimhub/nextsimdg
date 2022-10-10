/*!
 * @file ModelComponent_test.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "../src/include/ModelArrayRef.hpp"
#include "../src/include/ModelComponent.hpp"

#include <stdexcept>

namespace Nextsim {

// (Ab)use the exception mechanism to inform Catch that things are working correctly internally.
class HappyExcept : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

class Module1 : public ModelComponent {
public:
    Module1() { registerModule(); }
    std::string getName() const override { return "Module1"; }
    void setData(const ModelState::DataMap& st) override
    {
        throw(HappyExcept(std::string("setData for ") + getName()));
    }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }
    std::unordered_set<std::string> uFields() const override { return { "u1" }; }
    std::unordered_set<std::string> vFields() const override { return { "v1", "v2" }; }
    std::unordered_set<std::string> zFields() const override { return { "z1", "z2", "z3" }; }
};

TEST_CASE("Register a new module", "[ModelComponent]")
{
    Module1 m1;
    REQUIRE_THROWS_AS(ModelComponent::setAllModuleData(ModelState()), HappyExcept);

    std::unordered_set<std::string> uu;
    std::unordered_set<std::string> vv;
    std::unordered_set<std::string> zz;

    ModelComponent::getAllFieldNames(uu, vv, zz);
    REQUIRE(uu.size() == 1);
    REQUIRE(vv.size() == 2);
    REQUIRE(zz.size() == 3);

    ModelComponent::unregisterAllModules();
}

class ModuleSupplyAndWait : public ModelComponent {
public:
    ModuleSupplyAndWait()
        : hice(ModelArray::HField())
        , cice_ref(getProtectedArray())
    {
        registerModule();
        registerProtectedArray(ProtectedArray::H_ICE, &hice);
    }
    void setData(const ModelState::DataMap& ms) override { hice[0] = hiceData; }
    std::string getName() const override { return "SupplyAndWait"; }
    ModelState getState() const override
    {
        return {{
            { "hice", hice },
        }, {}};
    }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    const double hiceData = 1.2;
    double data() { return hice[0]; }
    double refData() { return cice_ref[0]; }

private:
    HField hice;
    ModelArrayRef<ProtectedArray::C_ICE, MARConstBackingStore> cice_ref;
};

class ModuleRequestAndSupply : public ModelComponent {
public:
    ModuleRequestAndSupply()
        : cice(ModelArray::HField())
        , hice_ref(getProtectedArray())
    {
        registerModule();
        registerProtectedArray(ProtectedArray::C_ICE, &cice);
    }
    void setData(const ModelState::DataMap& ms) override { cice[0] = ciceData; }
    std::string getName() const override { return "SupplyAndWait"; }
    ModelState getState() const override
    {
        return {{
            { "cice", cice },
        }, {}};
    }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    const double ciceData = 0.6;
    double data() { return cice[0]; }
    double refData() { return hice_ref[0]; }

private:
    HField cice;
    ModelArrayRef<ProtectedArray::H_ICE, MARConstBackingStore> hice_ref;
};

TEST_CASE("Test array registration", "[ModelComponent]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });
    ModuleSupplyAndWait saw;
    ModuleRequestAndSupply ras;

    REQUIRE(ras.data() == saw.refData());
    REQUIRE(saw.data() == ras.refData());
}

class ModuleSemiShared : public ModelComponent {
public:
    ModuleSemiShared()
        : qic(ModelArray::HField())
        , qio_ref(getSharedArray())
    {
        registerModule();
        registerSharedArray(SharedArray::Q_IC, &qic);
    }
    void setData(const ModelState::DataMap& ms) override { qic[0] = qicData; }
    std::string getName() const override { return "SemiShared"; }
    ModelState getState() const override
    {
        return {{
            { "qic", qic },
        }, {}};
    }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    const double qicData = 123;
    double data() { return qic[0]; }
    double refData() { return qio_ref[0]; }

private:
    HField qic;
    ModelArrayRef<SharedArray::Q_IO, MARBackingStore, RO> qio_ref;
};

class ModuleShared : public ModelComponent {
public:
    ModuleShared()
        : qio(ModelArray::HField())
        , qic_ref(getSharedArray())
    {
        registerModule();
        registerSharedArray(SharedArray::Q_IO, &qio);
    }
    void setData(const ModelState::DataMap& ms) override { qio[0]; }
    std::string getName() const override { return "Shared"; }
    ModelState getState() const override
    {
        return {{
            { "qio", qio },
        }, {}};
    }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    const double qioData = 234;
    const double qicData = 246;
    double data() { return qio[0]; }
    double& refData() { return qic_ref[0]; }
    void setRefData() { qic_ref[0] = qicData; }

private:
    HField qio;
    ModelArrayRef<SharedArray::Q_IC, MARBackingStore, RW> qic_ref;
};

TEST_CASE("Shared and semi-protected arrays", "[ModelComponent]")
{
    ModelArray::setDimensions(ModelArray::Type::H, { 1 });

    ModuleSemiShared semi;
    ModuleShared share;

    REQUIRE(share.data() == semi.refData());
    REQUIRE(semi.data() == share.refData());

    share.refData() = share.qicData;

    REQUIRE(semi.data() == share.qicData);
}

} /* namespace Nextsim */
