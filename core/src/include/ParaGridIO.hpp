/*!
 * @file ParaGridIO.hpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PARAGRIDIO_HPP
#define PARAGRIDIO_HPP

#include "include/ParametricGrid.hpp"

#include <map>
#include <ncFile.h>
#include <string>

namespace Nextsim {

/*!
 * A class to perform input and output for the ParametricGrid class. unlike the
 * other GridIO classes, this will hold non-restart files open to accumulate
 * data until closed using the close function.
 */
class ParaGridIO : public ParametricGrid::IParaGridIO {
public:
    ParaGridIO(ParametricGrid& grid)
        : IParaGridIO(grid)
    {
        if (dimCompMap.size() == 0)
            makeDimCompMap();
    }
    virtual ~ParaGridIO();

    /*!
     * Retrieves the ModelState from a restart file of the parametric_grid type.
     * @param filePath The file path containing the file to be read.
     */
    ModelState getModelState(const std::string& filePath) override;

    /*!
     * @brief Writes the ModelState to a given file location.
     *
     * @details Writes restart (isRestart true) or diagnostic (isRestart false)
     * output files from the provided model data and metadata. Restart files
     * contain a fixed set of fields without a time axis. Diagnostic files
     * output all provided fields and have a time axis to output data
     * throughout the model run.
     *
     * @params state The model state and configuration object.
     * @params metadata The model metadata (principally the initial file
     * creation model time).
     * @params filePath The path to output the data to. Identifies the output
     * stream for adding later time samples of the data fields.
     * @params isRestart Whether the file is a restart file or a diagnostic
     * output file.
     */
    void dumpModelState(const ModelState& state, const ModelMetadata& meta,
        const std::string& filePath, bool isRestart) override;

    /*!
     * Closes an open diagnostic file. Does nothing when provided with a
     * restart file name.
     *
     * @param filePath The path to the file to be closed.
     */
    void close(const std::string& filePath);

private:
    ParaGridIO() = delete;
    ParaGridIO(const ParaGridIO& other) = delete;
    ParaGridIO& operator=(const ParaGridIO& other) = delete;

    static const std::map<std::string, ModelArray::Type> dimensionKeys;

    static const std::map<ModelArray::Dimension, bool> isDG;
    static std::map<ModelArray::Dimension, ModelArray::Type> dimCompMap;

    // Ensures that static variables are created in the correct order.
    static void makeDimCompMap();

    std::map<std::string, netCDF::NcFile> openFiles;
    std::map<std::string, size_t> timeIndexByFile;
};

} /* namespace Nextsim */

#endif /* PARAGRIDIO_HPP */
