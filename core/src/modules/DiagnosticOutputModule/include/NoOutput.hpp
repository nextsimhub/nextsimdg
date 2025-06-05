/*!
 * @file NoOutput.hpp
 *
 * @date 23 Oct 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef NOOUTPUT_HPP
#define NOOUTPUT_HPP

#include "include/IDiagnosticOutput.hpp"

namespace Nextsim {

class NoOutput : public IDiagnosticOutput {
public:
    NoOutput() = default;

    void setFilenamePrefix(const std::string& filePrefix) override {};

    void outputState(const ModelState& state, const ModelMetadata& meta) override {};

    // ModelComponent functions
    std::string getName() const override { return "NoOutput"; }
};

} /* namespace Nextsim */

#endif /* NOOUTPUT_HPP */
