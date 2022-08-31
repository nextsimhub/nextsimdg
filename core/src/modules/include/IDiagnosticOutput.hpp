/*!
 * @file IDiagnosticOutput.hpp
 *
 * @date May 25, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDIAGNOSTICOUTPUT_HPP
#define IDIAGNOSTICOUTPUT_HPP

#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"

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
     * @param meta The model metadata for the the given state.
     */
    virtual void outputState(const ModelState& state, const ModelMetadata& meta) = 0;
};
}
#endif /* IDIAGNOSTICOUTPUT_HPP */
