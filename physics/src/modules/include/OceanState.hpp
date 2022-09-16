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
     * @brief Updates the ocean state at the start of the timestep.
     *
     * @details Sets the coupled-in SST/S (if applicable) and updates Tf for
     * all implementations.
     *
     * @param tst The object containing the timestep start and duration times.
     */
    virtual void updateBefore(const TimestepTime& tst);

    /*!
     * @brief Updates the ocean state at the end of the timestep.
     *
     * @details Calculates the nudging fluxes and updates the slab ocean.
     *
     * @param tst The object containing the timestep start and duration times.
     */
    virtual void updateAfter(const TimestepTime& tst);

protected:
    HField sst; // Slab ocean temperature, ˚C
    HField sss; // Slab ocean salinity, PSU
    HField mld; // Mixed layer (also slab ocean) depth, m

    HField tf; // Freezing point, ˚C
    HField cpml; // bulk heat capacity of the mixed layer, J m⁻² K⁻¹

    HField qdw; // Slab ocean heat flux, W m⁻²
    HField fdw; // Slab ocean water flux kg s ⁻¹ m⁻²

    IFreezingPoint* tfImpl;

    virtual void updateSpecial(const TimestepTime&) = 0;

    void updateFreezingPoint(const TimestepTime&);
private:
    void updateFreezingPointI(size_t i, const TimestepTime&);
};

} /* namespace Nextsim */

#endif /* OCEANSTATE_HPP */
