/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef TOPAZOCEAN_HPP
#define TOPAZOCEAN_HPP

#include "include/Configured.hpp"
#include "include/IFreezingPoint.hpp"
#include "include/IIceOceanHeatFlux.hpp"
#include "include/IOceanBoundary.hpp"
#include "include/SlabOcean.hpp"

namespace Nextsim {

/*!
 * A class to provided forcings from pre-processed forcings files based on ERA5
 * data.
 */
class TOPAZOcean : public IOceanBoundary, public Configured<TOPAZOcean> {
public:
    TOPAZOcean();
    ~TOPAZOcean() = default;

    enum {
        FILEPATH_KEY,
    };

    void setData(const ModelState::DataMap&) override;
    std::string getName() const override { return "TOPAZOcean"; }

    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    void configure() override;
    ConfigMap getConfiguration() const override;

    void updateBefore(const TimestepTime&) override;
    void updateAfter(const TimestepTime&) override;
    ModelState getStatePrognostic() const override;
    ModelState getStateDiagnostic() const override;

    void setFilePath(const std::string& filePathIn);

    const ModelArray getData(const std::string& nsName, const TimePoint& time) const;

    static void setDirectory(const std::string& dir);
    static const std::string& getDirectory();
    static const std::string addDirectory(const std::string& file);
private:
    // Location of the TOPAZ files
    static std::string& fileDirectory()
    {
        static std::string dir = ".";
        return dir;
    }
    // Updates the freezing point of an element
    void updateTf(size_t i, const TimestepTime& tst);
    // Since the configuration is global, it makes sense for the file path to
    // be static.
    static std::string filePath;

    HField sstExt;
    HField sssExt;

    SlabOcean slabOcean;

    std::unique_ptr<IIceOceanHeatFlux> pIOHeatFlux;
    std::unique_ptr<IFreezingPoint> pFreezingPoint;
};

} /* namespace Nextsim */

#endif /* TOPAZOCEAN_HPP_ */
