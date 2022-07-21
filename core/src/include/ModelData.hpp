/*!
 * @file ModelData.hpp
 *
 * @date Jul 20, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELDATA_HPP
#define MODELDATA_HPP

#include "include/AtmosphereOceanState.hpp"
#include "include/ModelMetadata.hpp"
#include "include/PrognosticData.hpp"
#include "include/Time.hpp"

namespace Nextsim {

class Model;

/*!
 * A class encapsulating all the data used when the model is running, as
 * opposed to setup and configuration data.
 */
class ModelData {
public:
    ModelData() = default;
    virtual ~ModelData() = default;

    /*!
     * @brief Updates the atmosphere and ocean states from the configured modules.
     * @param tst The start and duration of the current timestep
     */
    void updateAtmosOceanState(const TimestepTime& tst);
    /*!
     * @brief Steps the model forward by one timestep, changing the prognostic data.
     * @param tst The start and duration of the current timestep
     */
    void stepPrognosticData(const TimestepTime& tst);
    /*!
     * @brief Performs any additional export of data that is not diagnostic output.
     * @param tst The start and duration of the current timestep
     */
    void exportData(const TimestepTime& tst) const;
    //! Gets the current state of the model data fields.
    ModelState getModelState() const;
    //! Returns a constant reference to the metadata.
    const ModelMetadata& getMetadata() const { return metadata; }
private:
    PrognosticData pData;
    ModelMetadata metadata;
    AtmosphereOceanState aoState;

    friend Model;
};

} /* namespace Nextsim */

#endif /* MODELDATA_HPP */
