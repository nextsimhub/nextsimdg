/*!
 * @file ModelMetadata.hpp
 *
 * @date Jun 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELMETADATA_HPP
#define MODELMETADATA_HPP

#include "include/ConfigMap.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"

#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace Nextsim {

class CommonRestartMetadata;
/*!
 * A class to hold the metadata pertaining to the model as a whole, both
 * constant and time varying values. Especially values required for data file
 * output.
 */
class ModelMetadata {
public:
    /*!
     * @brief Sets the initial or current model time
     *
     * @param time TimePoint instance encoding the current time.
     */
    inline void setTime(const TimePoint& time) { m_time = time; }
    /*!
     * @brief Increments the model time metadata value.
     *
     * @param step Duration of the time increment to add.
     */
    inline void incrementTime(const Duration& step) { m_time += step; }
    //! Returns the current model time.
    inline const TimePoint& time() const { return m_time; }

    //! Returns the string description of the model grid structure.
    const std::string& structureName() const;

    /*!
     * @brief Sets the configuration metadata.
     *
     * @param config The configuration metadata
     */
    inline void setConfig(const ConfigMap& config) { m_config = config; }

    // The metadata writer should be a friend
    friend CommonRestartMetadata;

#ifdef USE_MPI
    void setMpiMetadata(MPI_Comm comm);

    MPI_Comm mpiComm;
    int mpiSize = 0;
    int mpiMyRank = -1;
    int localCornerX, localCornerY, localExtentX, localExtentY, globalExtentX, globalExtentY;
#endif

private:
    TimePoint m_time;
    ConfigMap m_config;

    // position coordinates on vertices
    ModelArray m_vertexCoords;
    // position coordinates of elements
    ModelArray m_coord1;
    ModelArray m_coord2;
    // Angle from model reference to grid north (+y for grids) TODO: what for meshes? N/A?
    ModelArray m_gridAzimuth;
};

} /* namespace Nextsim */

#endif /* MODELMETADATA_HPP */
