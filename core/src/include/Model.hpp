/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Kacper Kornet <kk562@cam.ac.uk>
 */

#ifndef MODEL_HPP
#define MODEL_HPP

#include "include/Logged.hpp"

#include "include/Configured.hpp"
#include "include/Iterator.hpp"
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

    void configureRestarts();
    void configureTime();
    void configure() override;

    enum {
        RESTARTFILE_KEY,
        STARTTIME_KEY,
        STOPTIME_KEY,
        RUNLENGTH_KEY,
        TIMESTEP_KEY,
        MISSINGVALUE_KEY,
#ifdef USE_MPI
        PARTITIONFILE_KEY,
#endif
        // Other Model configuration keys, not to be written to the restart file.
        RESTARTPERIOD_KEY,
        RESTARTOUTFILE_KEY,
#ifdef USE_OASIS
        WRITEOASISGRID_KEY,
#endif
    };

    ConfigMap getConfig() const { return ConfigMap(); };

    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    //! Run the model
    void run();

    void writeRestartFile();

    //! Sets the filename of the restart file that would currently be written out.
    void setFinalFilename(const std::string& finalFile);

    // Configuration option that holds the restart file name
    const static std::string restartOptionName;

private:
    Iterator iterator;
    DevStep modelStep; // Change the model step calculation here
    PrognosticData pData;
};

} /* namespace Nextsim */

#endif /* MODEL_HPP */
