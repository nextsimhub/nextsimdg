/*!
 * @file   ParaGridIO_Xios.hpp
 *
 * @date   27 Aug 2024
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 */
#ifdef USE_XIOS
#ifndef PARAGRIDIO_HPP
#define PARAGRIDIO_HPP

#include "ModelArray.hpp"
#include "StructureModule/include/ParametricGrid.hpp"

namespace Nextsim {

class ParaGridIO : public ParametricGrid::IParaGridIO {

public:
    ParaGridIO(ParametricGrid& grid)
        : IParaGridIO(grid)
    {
        // TODO: Implement this method
    }
    virtual ~ParaGridIO();

    // TODO: Align the API with ParaGridIO
    void read(const std::string fileId, ModelArray& modelarray);
    void write(const std::string fileId, ModelArray& modelarray);

    // TODO: Align the API with ParaGridIO
    /*!
     * Retrieves the ModelState from a restart file of the parametric_grid type.
     * @param filePath The file path containing the file to be read.
     */
#ifdef USE_MPI
    ModelState getModelState(const std::string& filePath, ModelMetadata& metadata) override;
#else
    ModelState getModelState(const std::string& filePath) override;
#endif

    // TODO: Align the API with ParaGridIO
    /*!
     * @brief Writes the ModelState to a given file location from the provided
     * model data and metadata.
     *
     * @params state The model state and configuration object.
     * @params metadata The model metadata (principally the initial file
     * creation model time).
     * @params filePath The path for the restart file.
     */
    void dumpModelState(
        const ModelState& state, const ModelMetadata& meta, const std::string& filePath) override;

    // TODO: Align the API with ParaGridIO
    /*!
     * @brief Reads forcings from a ParameticGrid flavoured file.
     *
     * @param forcings The names of the forcings required.
     * @param time The time for which to get the forcings.
     * @param filePath Path to the file to read.
     */
    ModelState readForcingTime(const std::set<std::string>& forcings, const TimePoint& time,
        const std::string& filePath) override;

    // TODO: Align the API with ParaGridIO
    /*!
     * @brief Writes diagnostic data to a file.
     *
     * @param state The state to write to the file.
     * @param time The time of the passed data.
     * @param filePath Path of the file to write to.
     */
    void writeDiagnosticTime(
        const ModelState& state, const ModelMetadata& meta, const std::string& filePath) override;
};
} /* namespace Nextsim */

#endif /* PARAGRIDIO_HPP */
#endif /* USE_XIOS */
