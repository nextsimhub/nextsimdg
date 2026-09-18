/*!
 *
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGOUTPUT_HPP
#define CONFIGOUTPUT_HPP

#include "include/IDiagnosticOutput.hpp"

#include "include/Configured.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"
#include "include/VectorRotator.hpp"
#include "include/gridNames.hpp"

#include <set>

namespace Nextsim {

/*!
 * An implementation of the diagnostic output that allows some configuration of
 * the file output period and frequency, as well as the fields the files contain.
 */
class ConfigOutput : public IDiagnosticOutput, public Configured<ConfigOutput> {
public:
    ConfigOutput();
    virtual ~ConfigOutput() = default;

    enum {
        PERIOD_KEY,
        START_KEY,
        SNAPSHOT_KEY,
        FIELDNAMES_KEY,
        FILENAME_KEY,
        FILEPERIOD_KEY,
        ORIENTATION_KEY,
    };

    // IDiagnosticOutput overrides
    void setFilenamePrefix(const std::string& filePrefix) override { m_filePrefix = filePrefix; }
    void setData(const TimePoint& modelStart) override;
    void outputState(const ModelState& diagState) override;

    // ModelComponent overrides
    inline std::string getName() const override { return "ConfigOutput"; };
    inline void setData(const ModelState::DataMap&) override {};

    // Configured overrides
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);
    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    void configure() override;

private:
    std::string m_filePrefix;
    Duration outputPeriod;
    bool firstOutput;
    bool everyTS;
    bool outputAllTheFields;
    TimePoint lastOutput;
    std::set<std::string> fieldsForOutput;
    std::string currentFileName;

    TimePoint lastFileChange;
    Duration fileChangePeriod;

    static const std::string all;
    static const std::string defaultLastOutput;

    std::map<std::string, std::string> reverseExternalNames;

    bool snapshots;
    bool resetState;

    std::vector<std::pair<std::string, std::string>> vectors
        = { { uName, vName }, { uWindName, vWindName }, { uOceanName, vOceanName } };
    std::unique_ptr<VectorRotator> rotator;
};

} /* namespace Nextsim */

#endif /* CONFIGOUTPUT_HPP */
