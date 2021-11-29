/*!
 * @file Model.hpp
 * @date 12 Aug 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef SRC_INCLUDE_MODEL_HPP
#define SRC_INCLUDE_MODEL_HPP

#include "Logged.hpp"

#include "Iterator.hpp"

namespace Nextsim {

class Model : public Logged {
public:
    Model(); // TODO add arguments to pass the desired
             // environment and configuration to the model
    ~Model(); // Finalize the model. Collect data and so on.
    void run();

private:
    Iterator iterator;
    Iterator::Iterant* iterant; // FIXME smart pointer

    bool deleteIterant;
};

} /* namespace Nextsim */

#endif /* SRC_INCLUDE_MODEL_HPP */
