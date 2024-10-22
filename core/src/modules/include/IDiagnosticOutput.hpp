/*!
 * @file IDiagnosticOutput.hpp
 *
 * @date 2 Jul 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDIAGNOSTICOUTPUT_HPP
#define IDIAGNOSTICOUTPUT_HPP

#include "include/ModelComponent.hpp"
#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"

#include <string>

namespace Nextsim {
class IDiagnosticOutput : public ModelComponent {
public:
    IDiagnosticOutput()
        : externalNames({
    /*
     * Using a pair of .ipp files to allow definitions of the externally visible
     * names of the fields to defined outside of an actual source file, even if the
     * definition file has a slightly odd format.
     */
#include "include/ProtectedArrayNames.ipp"
#include "include/SharedArrayNames.ipp"
          })
    {
    }
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
    virtual void outputState(const ModelMetadata& meta) = 0;

    // Define some of the ModelComponent class functions
    // No data to be set
    void setData(const ModelState::DataMap& state) { }

    /*!
     * Sets the model start time, which implementations may use to time output events.
     *
     * @brief modelStart the TimePoint of the start of the model run.
     */
    virtual void setModelStart(const TimePoint& ModelStart) { }

protected:
    const std::map<std::string, std::string> externalNames;
};
}
#endif /* IDIAGNOSTICOUTPUT_HPP */
