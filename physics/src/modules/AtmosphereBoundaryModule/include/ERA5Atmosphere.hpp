/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ERA5ATMOSPHERE_HPP
#define ERA5ATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"

#include "include/Configured.hpp"
#include "include/IFluxCalculation.hpp"

namespace Nextsim {

/*!
 * A class to provided forcings from pre-processed forcings files based on ERA5
 * data.
 */
class ERA5Atmosphere : public IAtmosphereBoundary, public Configured<ERA5Atmosphere> {
public:
    ERA5Atmosphere();
    ~ERA5Atmosphere() = default;

    enum {
        FILEPATH_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "ERA5Atmosphere"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;
    ConfigMap getConfiguration() const override;

    //! Calculates the fluxes from the given values
    void update(const TimestepTime&) override;

    void setFilePath(const std::string& filePathIn);

private:
    // Since the configuration is global, it makes sense for the file path to
    // be static.
    static std::string filePath;

    ModelArrayAccessor<Protected::T_AIR, RW> tairAccessor;
    ModelArrayAccessor<Protected::DEW_2M, RW> tdewAccessor;
    ModelArrayAccessor<Protected::P_AIR, RW> pairAccessor;
    ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;
    ModelArrayAccessor<Protected::LW_IN, RW> lw_inAccessor;
    ModelArrayAccessor<Protected::WIND_SPEED, RW> windAccessor;

    IFluxCalculation* fluxImpl;
};

} /* namespace Nextsim */

#endif /* ERA5ATMOSPHERE_HPP_ */
