/*!
 * @file SimpleOutput.hpp
 *
 * @date 24 Sep 2024
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

    void outputState(const ModelMetadata& meta) override;

    // ModelComponent functions
    std::string getName() const override { return "SimpleOutput"; }
    // No configuration in getState
    ModelState getState() const override { return ModelState(); }
    ModelState getState(const OutputLevel&) const override { return getState(); }

private:
    std::string m_filePrefix;
};

} /* namespace Nextsim */

#endif /* SIMPLEOUTPUT_HPP */
