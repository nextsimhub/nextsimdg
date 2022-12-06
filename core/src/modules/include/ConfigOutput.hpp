/*!
 * @file ConfigOutput.hpp
 *
 * @date Aug 22, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CONFIGOUTPUT_HPP
#define CONFIGOUTPUT_HPP

#include "include/IDiagnosticOutput.hpp"

#include "include/Configured.hpp"
#include "include/Time.hpp"

#include <set>

namespace Nextsim {

/*!
 * An implementation of the diagnostic output that allows some configuration of
 * the file output period and frequency, as well as the fields the files contain.
 */
class ConfigOutput : public IDiagnosticOutput, public Configured<ConfigOutput> {
public:
    ConfigOutput() = default;
    virtual ~ConfigOutput() = default;

    enum {
        PERIOD_KEY,
        START_KEY,
        FIELDNAMES_KEY,
    };

    void setFilenamePrefix(const std::string& filePrefix) override { m_filePrefix = filePrefix; }

    void setOutputDirectory(const std::filesystem::path outputDirectory) override {
        m_outputDirectory = outputDirectory;
    }

    void outputState(const ModelState& state, const ModelMetadata& meta) override;

    void configure() override;

private:
    std::string m_filePrefix;
    std::filesystem::path m_outputDirectory;
    Duration outputPeriod;
    bool firstOutput = true;
    bool everyTS = false;
    bool outputAllTheFields = false;
    TimePoint lastOutput;
    std::set<std::string> fieldsForOutput;

    static const std::string all;
};

} /* namespace Nextsim */

#endif /* CONFIGOUTPUT_HPP */
