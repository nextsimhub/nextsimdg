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
    ModelState getState(Logged::level& lvl) const override { return getState(); }
};

TEST_CASE("Register a new module", "[ModelModule]")
{
    Module1 m1;
    REQUIRE_THROWS_AS(ModelModule::setAllModuleData(ModelState()), HappyExcept);
    ModelModule::unregisterAllModules();
}

} /* namespace Nextsim */
