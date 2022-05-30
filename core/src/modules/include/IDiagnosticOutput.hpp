/*!
 * @file IDiagnosticOutput.hpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDIAGNOSTICOUTPUT_HPP
#define IDIAGNOSTICOUTPUT_HPP

#include "include/ModelState.hpp"
#include "include/Time.hpp"

#include <string>

namespace Nextsim {
class IDiagnosticOutput {
public:
    IDiagnosticOutput() = default;
    virtual ~IDiagnosticOutput() = default;

    /*!
     * @brief Sets the output file name.
     *
     * @param fileName The file name to be set.
     */
    virtual void setFilenamePrefix(const std::string& filePrefix) = 0;

    /*!
     * @brief Outputs the passed ModelState.
     *
     * @param state The model state to be written out.
     * @param tst The time data of the current timestep.
     */
    virtual void outputState(const ModelState& state, const TimestepTime& tst) const = 0;
};
}
#endif /* IDIAGNOSTICOUTPUT_HPP */
