/*!
 * @file ParaGridIO_Xios.hpp
 *
 * @date 2 August, 2024
 * @author Joe Wallwork <jw2423@cam.ac.uk>
 */

#ifndef PARAGRIDIO_XIOS_HPP
#define PARAGRIDIO_XIOS_HPP

#include "StructureModule/include/ParametricGrid.hpp"

namespace Nextsim {

/*!
 * A class to perform input and output for the ParametricGrid class using XIOS
 * asynchronous input-output tool, replicating the functionality of ParaGridIO.
 * This class also replicates its behaviour in that it will hold non-restart
 * files open to accumulate data until the close function.
 */
class ParaGridIO : public ParametricGrid::IParaGridIO {
public:
    ParaGridIO(ParametricGrid& grid)
        : IParaGridIO(grid)
    {
        // TODO: Implement this method
    }
    virtual ~ParaGridIO();

    /*!
     * Retrieves the ModelState from a restart file of the ParametricGrid type.
     *
     * @param filePath The file path containing the file to be read.
     */
    // TODO: Implement this method
    ModelState getModelState(const std::string& filePath) override;

    /*!
     * @brief Writes the ModelState to a given file location from the provided
     * model data and metadata.
     *
     * @param state The model state and configuration object.
     * @param metadata The model metadata (principally the initial file
     * creation model time).
     * @param filePath The path for the restart file.
     */
    // TODO: Implement this method
    void dumpModelState(
        const ModelState& state, const ModelMetadata& meta, const std::string& filePath) override;

    /*!
     * @brief Reads forcings from a ParameticGrid flavoured file.
     *
     * @param forcings The names of the forcings required.
     * @param time The time for which to get the forcings.
     * @param filePath Path to the file to read.
     */
    ModelState readForcingTime(const std::set<std::string>& forcings, const TimePoint& time,
        const std::string& filePath) override
    {
        return readForcingTimeStatic(forcings, time, filePath);
    }

    /*!
     * @brief Writes diagnostic data to a file.
     *
     * @param state The state to write to the file.
     * @param time The time of the passed data.
     * @param filePath Path of the file to write to.
     */
    // TODO: Implement this method
    void writeDiagnosticTime(
        const ModelState& state, const ModelMetadata& meta, const std::string& filePath) override;

    /*!
     * Closes an open diagnostic file. Does nothing when provided with a
     * restart file name.
     *
     * @param filePath The path to the file to be closed.
     */
    // TODO: Implement this method
    static void close(const std::string& filePath);

    // TODO: Implement this method
    static ModelState readForcingTimeStatic(
        const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath);

private:
    /* TODO: Implement private methods? */
};

} /* namespace Nextsim */

#endif /* PARAGRIDIO_XIOS_HPP */
