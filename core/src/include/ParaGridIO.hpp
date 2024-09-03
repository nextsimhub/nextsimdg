/*!
 * @file ParaGridIO.hpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PARAGRIDIO_HPP
#define PARAGRIDIO_HPP

#include "StructureModule/include/ParametricGrid.hpp"

#ifdef USE_MPI
#include "ParallelNetcdfFile.hpp"
#endif
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
#ifdef USE_MPI
    typedef NetCDFFileType netCDF::NcFilePar;
#else
    typedef netCDF::NcFile NetCDFFileType ;
#endif

    ParaGridIO(ParametricGrid& grid)
        : IParaGridIO(grid)
        , openFiles(getOpenFiles())
        , timeIndexByFile(getFileTimeIndices())
    {
        if (dimCompMap.size() == 0)
            makeDimCompMap();
    }
    virtual ~ParaGridIO();

    /*!
     * Retrieves the ModelState from a restart file of the parametric_grid type.
     * @param filePath The file path containing the file to be read.
     */
#ifdef USE_MPI
    ModelState getModelState(const std::string& filePath, ModelMetadata& metadata) override;
#else
    ModelState getModelState(const std::string& filePath) override;
#endif

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
    void writeDiagnosticTime(
        const ModelState& state, const ModelMetadata& meta, const std::string& filePath) override;

    /*!
     * Closes an open diagnostic file. Does nothing when provided with a
     * restart file name.
     *
     * @param filePath The path to the file to be closed.
     */
    static void close(const std::string& filePath);

    static ModelState readForcingTimeStatic(
        const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath);

private:
    typedef std::map<std::string, NetCDFFileType> FileMap;
    typedef std::map<std::string, size_t> IndexMap;

    ParaGridIO() = delete;
    ParaGridIO(const ParaGridIO& other) = delete;
    ParaGridIO& operator=(const ParaGridIO& other) = delete;

    static const std::map<std::string, ModelArray::Type> dimensionKeys;

    static const std::map<ModelArray::Dimension, bool> isDG;
    static std::map<ModelArray::Dimension, ModelArray::Type> dimCompMap;

    // Ensures that static variables are created in the correct order.
    static void makeDimCompMap();

    // Closes all still-open NetCDF files
    static void closeAllFiles();

    // Existing or open files are a property of the computer outside the individual
    // class instance, so they are singletons.
    FileMap& openFiles;
    inline static FileMap& getOpenFiles()
    {
        static FileMap fm;
        return fm;
    }

    IndexMap& timeIndexByFile;
    inline static IndexMap& getFileTimeIndices()
    {
        static IndexMap tm;
        return tm;
    }
};

} /* namespace Nextsim */

#endif /* PARAGRIDIO_HPP */
