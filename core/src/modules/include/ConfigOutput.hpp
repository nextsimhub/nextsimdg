/*!
 * @file ConfigOutput.hpp
 *
 * @date 7 Sep 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGOUTPUT_HPP
#define CONFIGOUTPUT_HPP

#include "include/IDiagnosticOutput.hpp"

#include "include/Configured.hpp"
#include "include/ModelComponent.hpp"
#include "include/Time.hpp"

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
        FIELDNAMES_KEY,
        FILENAME_KEY,
        FILEPERIOD_KEY,
    };

    // IDiagnosticOutput overrides
    void setFilenamePrefix(const std::string& filePrefix) override { m_filePrefix = filePrefix; }

    void outputState(const ModelMetadata& meta) override;

    // ModelComponent overrides
    inline std::string getName() const override { return "ConfigOutput"; };
    inline void setData(const ModelState::DataMap&) override {};
    inline ModelState getState() const override { return ModelState(); };
    inline ModelState getState(const OutputLevel&) const override { return ModelState(); };

    // Configured overrides
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);
    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    void configure() override;
    ModelState getStateRecursive(const OutputSpec& os) const override;

private:
    std::string m_filePrefix;
    Duration outputPeriod;
    bool firstOutput = true;
    bool everyTS = false;
    bool outputAllTheFields = false;
    TimePoint lastOutput;
    std::set<std::string> fieldsForOutput;
    std::string currentFileName;
    std::set<std::string> internalFieldsForOutput;

    TimePoint lastFileChange;
    Duration fileChangePeriod;

    TimePoint lastFileChange;
    Duration fileChangePeriod;

    static const std::string all;
    static const std::string defaultLastOutput;

    std::map<std::string, std::string> reverseExternalNames;

};

} /* namespace Nextsim */

#endif /* CONFIGOUTPUT_HPP */
