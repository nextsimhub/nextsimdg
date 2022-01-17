/*!
 * @file Model.hpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_MODEL_HPP
#define SRC_INCLUDE_MODEL_HPP

#include "include/Logged.hpp"

#include "include/Configured.hpp"
#include "include/Iterator.hpp"
#include "include/IStructure.hpp"

#include <string>
#include "DevStep.hpp"

namespace Nextsim {

//! A class that encapsulates the whole of the model
class Model : public Logged, public Configured<Model> {
public:
    Model(); // TODO add arguments to pass the desired
             // environment and configuration to the model
    ~Model(); // Finalize the model. Collect data and so on.

    void configure() override;
    enum {
        RESTARTFILE_KEY,
        STARTTIME_KEY,
        STOPTIME_KEY,
        RUNLENGTH_KEY,
        TIMESTEP_KEY,
    };

    //! Run the model
    void run();

private:
    Iterator iterator;
    DevStep iterant; // Change the model master iterant here

    std::string restartFileName;

    IStructure* dataStructure;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_MODEL_HPP */
