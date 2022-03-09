/*!
 * @file ModelModule_test.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/ModelModule.hpp"

#include <stdexcept>

namespace Nextsim {

// (Ab)use the exception mechanism to inform Catch that things are working correctly internally.
class HappyExcept : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

class Module1 : public ModelModule {
public:
    Module1() { registerModule(); }
    std::string getName() const override { return "Module1"; }
    void setData(const ModelState& st) override
    {
        throw(HappyExcept(std::string("setData for ") + getName()));
    }
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }
    std::set<std::string> uFields() const override { return { "u1" }; }
    std::set<std::string> vFields() const override { return { "v1", "v2" }; }
    std::set<std::string> zFields() const override { return { "z1", "z2", "z3" }; }
};

TEST_CASE("Register a new module", "[ModelModule]")
{
    Module1 m1;
    REQUIRE_THROWS_AS(ModelModule::setAllModuleData(ModelState()), HappyExcept);

    std::set<std::string> uu;
    std::set<std::string> vv;
    std::set<std::string> zz;

    ModelModule::getAllFieldNames(uu, vv, zz);
    REQUIRE(uu.size() == 1);
    REQUIRE(vv.size() == 2);
    REQUIRE(zz.size() == 3);

    ModelModule::unregisterAllModules();
}

class ModuleSupplyAndWait : public ModelModule {
public:
    ModuleSupplyAndWait()
        : hice(ModelArray::HField("hice"))
        , p_cice(nullptr)
    {
        registerModule();
        registerProtectedArray(ProtectedArray::H_ICE, &hice);
        requestProtectedArray(ProtectedArray::C_ICE, &p_cice);
    }
    void setData(const ModelState& ms) override { }
    std::string getName() const override { return "SupplyAndWait"; }
    ModelState getState() const override
    {
        return {
            { "hice", hice },
        };
    }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    bool checkNotNull() { return p_cice; }

private:
    HField hice;
    const HField* p_cice;
};

class ModuleRequestAndSupply : public ModelModule {
public:
    ModuleRequestAndSupply()
        : cice(ModelArray::HField("cice"))
        , p_hice(nullptr)
    {
        registerModule();
        registerProtectedArray(ProtectedArray::C_ICE, &cice);
        requestProtectedArray(ProtectedArray::H_ICE, &p_hice);
    }
    void setData(const ModelState& ms) override { }
    std::string getName() const override { return "SupplyAndWait"; }
    ModelState getState() const override
    {
        return {
            { "cice", cice },
        };
    }
    ModelState getState(const OutputLevel& lvl) const override { return getState(); }

    bool checkNotNull() { return p_hice; }

private:
    HField cice;
    const HField* p_hice;
};

TEST_CASE("Test array registration", "[ModelModule]")
{
    ModuleSupplyAndWait saw;
    ModuleRequestAndSupply ras;

    REQUIRE(ras.checkNotNull());
    REQUIRE(saw.checkNotNull());
}

} /* namespace Nextsim */
