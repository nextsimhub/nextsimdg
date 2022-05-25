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

    void setFilename(const std::string& fileName) override { m_fileName = fileName; }

    void outputState(const ModelState& state, const TimestepTime& tst) const override;

private:
    std::string m_fileName;
};

} /* namespace Nextsim */

#endif /* SIMPLEOUTPUT_HPP */
