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

    void outputState(const ModelState& state, const ModelMetadata& meta) const override;

    void configure() override;

private:
    std::string m_filePrefix;
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
