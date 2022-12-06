/*!
 * @file IModelStep.hpp
 *
 * @date Jan 18, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IMODELSTEP_HPP
#define IMODELSTEP_HPP

#include "include/Iterator.hpp"
#include "include/ModelMetadata.hpp"
#include "include/PrognosticData.hpp"

namespace Nextsim {

//! An interface class extending Iterator::Iterant with some file name handling.
class IModelStep : public Iterator::Iterant {
public:
    IModelStep() = default;
    virtual ~IModelStep() = default;
    /*!
     * @brief Sets the path to the initial file for later reference.
     *
     * @param filePath The path to the location of the file used to initialize
     *        the model.
     */
    void setInitFile(const std::string& filePath) { initialRestartFilePath = filePath; };

    /*!
     * @brief Sets the path to the output directory for later reference.
     *
     * @param filePath The path to the location in which diagnostic files and model outputs
     *        should be stored.
     */
    //TODO: See if we can do without this member variable before submission
    void setOutputDir(const std::filesystem::path outputDir) { outputDirectory = outputDir; };

    /*!
     * @brief Writes a restart file containing the current model state to the
     * given file location
     *
     * @param filePath The path to the location to write the file.
     */
    virtual void writeRestartFile(const std::string& filePath) = 0;

    /*!
     * @brief Sets the data object that will be used within the timesteps.
     *
     * @param data The PrognosticData holding the model data.
     */
    virtual void setData(PrognosticData& data) = 0;

    /*!
     * @brief Sets the metadata object that will be used within the timesteps.
     *
     * @details The ModelMetadata object will contain much constant
     * information, primarily metadata about the model set up, but also
     * time-dependent metadata, including the current time for which the model
     * data is valid.
     *
     * @param meta A ModelMetadata reference to the metadata instance
     */
    virtual void setMetadata(ModelMetadata& meta) = 0;

    // Member functions inherited from Iterant
    virtual void init() = 0;
    virtual void start(const TimePoint& startTime) = 0;
    virtual void iterate(const TimestepTime& dt) = 0;
    virtual void stop(const TimePoint& stopTime) = 0;

protected:
    std::string initialRestartFilePath;
    std::filesystem::path outputDirectory;
};

} /* namespace Nextsim */

#endif /* IMODELSTEP_HPP */
