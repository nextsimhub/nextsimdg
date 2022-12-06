/*!
 * @file SimpleOutput.hpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SIMPLEOUTPUT_HPP
#define SIMPLEOUTPUT_HPP

#include "include/IDiagnosticOutput.hpp"

namespace Nextsim {

class SimpleOutput : public IDiagnosticOutput {
public:
    SimpleOutput() = default;

    void setFilenamePrefix(const std::string& filePrefix) override { m_filePrefix = filePrefix; }
    void setOutputDirectory(const std::filesystem::path outputDirectory) override {
        m_outputDirectory = outputDirectory;
    }

    void outputState(const ModelState& state, const ModelMetadata& meta) override;

private:
    std::string m_filePrefix;
    std::filesystem::path m_outputDirectory;
};

} /* namespace Nextsim */

#endif /* SIMPLEOUTPUT_HPP */
