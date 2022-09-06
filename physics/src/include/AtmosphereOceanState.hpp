/*!
 * @file AtmosphereOceanState.hpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ATMOSPHEREOCEANSTATE_HPP
#define ATMOSPHEREOCEANSTATE_HPP

#include "include/AtmosphereState.hpp"
#include "include/Configured.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArrayRef.hpp"
#include "include/ModelComponent.hpp"

#include <memory>

namespace Nextsim {

class AtmosphereOceanState : public ModelComponent, public Configured<AtmosphereOceanState> {
public:
    AtmosphereOceanState();

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

    std::string getName() const override;
    std::unordered_set<std::string> hFields() const override;

    void update(const TimestepTime&);

    void configure() override;

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

protected:
    std::unique_ptr<AtmosphereState> atmosStateImpl;
};

} /* namespace Nextsim */

#endif /* ATMOSPHEREOCEANSTATE_HPP */
