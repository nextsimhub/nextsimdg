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
        DIRPATH_KEY,
        CHECKOVERRIDE_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "ERA5Atmosphere"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;
    ConfigMap getConfiguration() const override;

    //! Calculates the fluxes from the given values
    void update(const TimestepTime&) override;

    void setFilePath(const std::string& filePathIn);

    static void setDirectory(const std::string& dir);
    static const std::string& getDirectory();
    static const std::string addDirectory(const std::string& file);

private:
    // Since the configuration is global, it makes sense for the file path to
    // be static.
    static std::string filePath;
    // Location of the ERA5 files
    static std::string& fileDirectory()
    {
        static std::string dir = ".";
        return dir;
    }
    static bool& checkOverride()
    {
        static bool ignore = false;
        return ignore;
    }

    const ModelArray getData(const std::string& nsName, const TimePoint& time) const;

    HField tair;
    HField tdew;
    HField pair;
    HField sw_in;
    HField lw_in;
    HField wind;

    HField modelLon;
    HField modelLat;

    IFluxCalculation* fluxImpl;
};

} /* namespace Nextsim */

#endif /* ERA5ATMOSPHERE_HPP_ */
