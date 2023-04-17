/*!
 * @file Model.hpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODEL_HPP
#define MODEL_HPP

#include "include/Logged.hpp"

#include "include/Configured.hpp"
#include "include/Iterator.hpp"
#include "include/ModelConfig.hpp"
#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"
#include "include/PrognosticData.hpp"

#include "DevStep.hpp"
#include <string>

namespace Nextsim {

//! A class that encapsulates the whole of the model
class Model : public Configured<Model> {
public:
    Model(); // TODO add arguments to pass the desired
             // environment and configuration to the model
    ~Model(); // Finalize the model. Collect data and so on.

    void configure() override;

    ConfigMap getConfig() const;

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    //! Run the model
    void run();

    void writeRestartFile();

    //! Sets the filename of the restart file that would currently be written out.
    void setFinalFilename(const std::string& finalFile);

    //! Gets the model metadata instance
    ModelMetadata& metadata();

    // Configuration option that holds the restart file name
    const static std::string restartOptionName;

    enum {
        RESTARTFILE_KEY,
        STARTTIME_KEY = ModelConfig::STARTTIME_KEY,
        STOPTIME_KEY = ModelConfig::STOPTIME_KEY,
        RUNLENGTH_KEY = ModelConfig::RUNLENGTH_KEY,
        TIMESTEP_KEY = ModelConfig::TIMESTEP_KEY,
        RESTARTPERIOD_KEY = ModelConfig::RESTARTPERIOD_KEY,
    };

private:
    Iterator iterator;
    DevStep modelStep; // Change the model step calculation here
    PrognosticData pData;
    ModelMetadata m_etadata;

    std::string initialFileName;
    std::string finalFileName;
};

} /* namespace Nextsim */

#endif /* MODEL_HPP */
