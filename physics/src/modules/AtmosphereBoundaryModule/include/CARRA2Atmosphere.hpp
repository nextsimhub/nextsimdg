/*!
 * @author  Einar Ólason <einar.olason@nersc.no>
 */

#ifndef CARRA2ATMOSPHERE_HPP
#define CARRA2ATMOSPHERE_HPP

#include "include/IAtmosphereBoundary.hpp"

#include "include/Configured.hpp"
#include "include/IFluxCalculation.hpp"
#include "include/IForcingInput.hpp"

namespace Nextsim {

/*!
 * A class to provided forcings from pre-processed forcings files based on CARRA2
 * data.
 */
class CARRA2Atmosphere : public IAtmosphereBoundary, public Configured<CARRA2Atmosphere> {
public:
    CARRA2Atmosphere();

    enum {
        FILEPATH_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "CARRA2Atmosphere"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;

    //! Calculates the fluxes from the given values
    void update(const TimestepTime&) override;

    static void setFilePath(const std::string& filePathIn);

protected:
    ConfigMap getConfiguration() const override;

private:
    // Since the configuration is global, it makes sense for the file path to
    // be static.
    static std::string filePath;

    // A list of variable names for the CARRA2 inputs
    const std::string lwInName = "strd";
    const std::string pAirName = "msl";
    const std::string rainName = "tirf";
    const std::string snowName = "titspf";
    const std::string sphumName = "2sh";
    const std::string swInName = "ssrd";
    const std::string tAirName = "2t";
    const std::string uName = "10u";
    const std::string vName = "10v";

    const std::string ncLonName = "lon";
    const std::string ncLatName = "lat";
    const std::string ncTimeName = "time";

    const std::set<std::string> forcings
        = { sphumName, lwInName, pAirName, rainName, snowName, swInName, tAirName, uName, vName };
    const std::set<std::pair<std::string, std::string>> vectors = { { uName, vName } };

    std::unique_ptr<IForcingInput> forcingState;

    ModelArrayAccessor<Protected::T_AIR, RW> tairAccessor;
    ModelArrayAccessor<Protected::SPHUM, RW> spHumAccessor;
    ModelArrayAccessor<Protected::P_AIR, RW> pairAccessor;
    ModelArrayAccessor<Protected::SW_IN, RW> sw_inAccessor;
    ModelArrayAccessor<Protected::LW_IN, RW> lw_inAccessor;
    ModelArrayAccessor<Protected::WIND_SPEED, RW> windAccessor;

    std::unique_ptr<IFluxCalculation> fluxImpl;
};

} /* namespace Nextsim */
#endif // CARRA2ATMOSPHERE_HPP
