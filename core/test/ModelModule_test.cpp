/*!
 * @file ModelModule_test.cpp
 *
 * @date Feb 28, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "include/ModelModule.hpp"

namespace Nextsim {

class Module1 : public ModelModule {
public:
    Module1()
    {
        registerModule();
    }
    std::string getName() const override
    {
        return "Module1";
    }
    void setData(const ModelState& st) override
    {
        INFO("setData for " << getName());
    }
    ModelState getState() const override
    {
        return ModelState();
    }
    ModelState getState(Logged::level& lvl) const override
    {
        return getState();
    }
};

TEST_CASE("Register a new module", "[ModelModule]")
{
    Module1 m1;
    ModelModule::setAllModuleData(ModelState());
}

} /* namespace Nextsim */
