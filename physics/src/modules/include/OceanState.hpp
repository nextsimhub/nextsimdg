/*!
 * @file OceanState.hpp
 *
 * @date May 9, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef OCEANSTATE_HPP
#define OCEANSTATE_HPP

#include "include/Configured.hpp"
#include "include/IFreezingPointModule.hpp"
#include "include/ModelComponent.hpp"

namespace Nextsim {

//! A class providing an interface to the ocean climatologies or coupled model.
class OceanState : public ModelComponent, public Configured<OceanState> {
public:
    OceanState();
    virtual ~OceanState() = default;

    void setData(const ModelState::DataMap&) override;
    ModelState getState() const override;
    ModelState getState(const OutputLevel&) const override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

    std::string getName() const override;
    std::unordered_set<std::string> hFields() const override;

    void configure() override;

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);
    /*!
     * @brief Updates the ocean state.
     *
     * @details Performs any common calculations, then any implementation
     * specific updates.
     *
     * @param tStep The object containing the timestep start and duration times.
     */
    void update(const TimestepTime&);

    //! Returns whether the implementation uses the slab ocean (defaults to true)
    virtual bool usesSlabOcean() { return b_usesSlabOcean; }

protected:
    HField sst;
    HField sss;
    HField mld;

    HField tf;
    HField cpml;

    IFreezingPoint* tfImpl;

    bool b_usesSlabOcean;
    virtual void updateSpecial(const TimestepTime&) = 0;

private:
    void updateFreezingPoint(size_t i, const TimestepTime&);
};

} /* namespace Nextsim */

#endif /* OCEANSTATE_HPP */
