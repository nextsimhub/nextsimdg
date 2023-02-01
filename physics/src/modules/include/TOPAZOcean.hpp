/*!
 * @file TOPAZOcean.hpp
 *
 * @date Nov 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef TOPAZOCEAN_HPP
#define TOPAZOCEAN_HPP

#include "IOceanBoundary.hpp"

#include "include/Configured.hpp"
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

    void updateBefore(const TimestepTime&) override;
    void updateAfter(const TimestepTime&) override;

    void setFilePath(const std::string& filePathIn);

private:
    // Updates the freezing point of an element
    void updateTf(size_t i, const TimestepTime& tst);
    // Since the configuration is global, it makes sense for the file path to
    // be static.
    static std::string filePath;

    SlabOcean slabOcean;
};

} /* namespace Nextsim */

#endif /* TOPAZOCEAN_HPP_ */
