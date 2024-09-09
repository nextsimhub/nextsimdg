/*!
 * @file ModelMetadata.hpp
 *
 * @date 09 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MODELMETADATA_HPP
#define MODELMETADATA_HPP

#include "include/ConfigMap.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"

#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_OASIS
#include <oasis_c.h>
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
#ifdef USE_MPI
    /*!
     * @brief Construct a ModelMedata based on file with decompostion data
     *
     * @param partitionFile The name of file with decompostion data
     * @param comm MPI communicator
     */
    ModelMetadata(std::string partitionFile, MPI_Comm comm);

    // We need to force default constructor also
    ModelMetadata() = default;
#endif

#ifdef USE_OASIS
    int OASISPartitionId;
#endif

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

    /*!
     * @brief Extracts and sets the coordinate metadata from the given ModelState.
     *
     * @param state The given ModelState.
     */
    const ModelState& extractCoordinates(const ModelState& state);

    /*!
     * @brief Adds the coordinate metadata to the given ModelState.
     *
     * @param state The given ModelState.
     */
    ModelState& affixCoordinates(ModelState& state) const;

#ifdef USE_MPI
    void setMpiMetadata(MPI_Comm comm);
    /*!
     * @brief Extracts and sets MPI partition metadata from partition file
     *
     * @param file with partition metadata
     */
    void getPartitionMetadata(std::string partitionFile);

    MPI_Comm mpiComm;
    int mpiSize = 0;
    int mpiMyRank = -1;
    int localCornerX, localCornerY, localExtentX, localExtentY, globalExtentX, globalExtentY;
#endif

#ifdef USE_OASIS
    void initOasis(const bool writeOasisGrid)
    {
        // Set the partitioning
        /* From the manual: "[ig_paral is a] vector of integers describing the local grid partition
         * in the global index space; has a different expression depending on the type of the
         * partition; in OASIS3-MCT, 5 types of partition are supported: Serial (no partition),
         * Apple, Box, Orange, and Points" - it looks like we should use "Box", so partInfo[0] = 2
         * (aka. ig_paral).
         *
         * #Box partition#
         * Each partition is a rectangular region of the global domain, described by the global
         * offset of its upper left corner, and its local extents in the X and Y dimensions. The
         * global extent in the X dimension must also be given. In this case, we have ig_paral(1:5):
         *  - ig_paral(1) = 2 (indicates a Box partition)
         *  - ig_paral(2) = the upper left corner global offset
         *  - ig paral(3) = the local extent in x
         *  - ig_paral(4) = the local extent in y
         *  - ig_paral(5) = the global extent in x.
         *
         * metdatata contains: localCornerX, localCornerY, localExtentX, localExtentY,
         * globalExtentX, globalExtentY;
         */
        // TODO: The contents of metadata is not certain!
        const int offset = globalExtentX * localCornerY + localCornerX;
        const std::vector<int> partInfo
            = { OASIS_Box, offset, localExtentX, localExtentY, globalExtentX };

        const int globalSize = globalExtentX * globalExtentY;
        const std::string compName = "nextsim"; // Not useful for any setups we have in mind
        OASIS_CHECK_ERR(oasis_c_def_partition(
            &OASISPartitionId, OASIS_Box_Params, &partInfo[0], globalSize, compName.c_str()));

        // TODO: Writing out grid information should be possible, but optional
        if (writeOasisGrid) { }
    }
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
    // Are the coordinates Cartesian? x & y versus longitude and latitude
    bool isCartesian;
    // Are the more complex coordinates stored?
    bool hasParameters;
#ifdef USE_MPI
    const std::string bboxName = "bounding_boxes";
#endif
};

} /* namespace Nextsim */

#endif /* MODELMETADATA_HPP */
