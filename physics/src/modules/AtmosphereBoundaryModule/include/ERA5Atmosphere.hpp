/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ERA5ATMOSPHERE_HPP
#define ERA5ATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"

#include "include/Configured.hpp"
#include "include/IFluxCalculation.hpp"
#include "include/ParaGridInputs.hpp"

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

    static void setFilePath(const std::string& filePathIn);

private:
    // Since the configuration is global, it makes sense for the file path to
    // be static.
    static std::string filePath;

    // A list of variable names for the ERA5 inputs
    const std::string tAirName = "t2m";
    const std::string dew2mName = "d2m";
    const std::string pAirName = "msl";
    const std::string swInName = "msdwswrf";
    const std::string lwInName = "msdwlwrf";
    const std::string snowName = "msr";
    const std::string rainName = "mtpr";
    const std::string uName = "u10";
    const std::string vName = "v10";

    const std::string ncLonName = "longitude";
    const std::string ncLatName = "latitude";
    const std::string ncTimeName = "time";

    const std::set<std::string> forcings
        = { tAirName, dew2mName, pAirName, swInName, lwInName, uName, vName };
    const std::set<std::pair<std::string, std::string>> vectors = { { uName, vName } };

    ParaGridInputs forcingState;

    ModelArrayAccessor<Protected::T_AIR, RW> tairAccessor;
    ModelArrayAccessor<Protected::DEW_2M, RW> tdewAccessor;
    ModelArrayAccessor<Protected::P_AIR, RW> pairAccessor;
    ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;
    ModelArrayAccessor<Protected::LW_IN, RW> lw_inAccessor;
    ModelArrayAccessor<Protected::WIND_SPEED, RW> windAccessor;

    std::unique_ptr<IFluxCalculation> fluxImpl;
};

} /* namespace Nextsim */

#endif /* ERA5ATMOSPHERE_HPP_ */
