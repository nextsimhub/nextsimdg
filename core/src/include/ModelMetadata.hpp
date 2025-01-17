/*!
 * @file ModelMetadata.hpp
 *
 * @date 17 Jan 2025
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifndef MODELMETADATA_HPP
#define MODELMETADATA_HPP

#include "include/ConfigMap.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelState.hpp"
#include "include/Time.hpp"
#include <ncFile.h>
#include <string>
#include <vector>

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

    /*!
     * @brief Sets the initial or current model time
     *
     * @param time TimePoint instance encoding the current time.
     */
    void setTime(const TimePoint& time);
    /*!
     * @brief Increments the model time metadata value.
     *
     * @param step Duration of the time increment to add.
     */
    void incrementTime(const Duration& step);
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

    enum Edge { BOTTOM, RIGHT, TOP, LEFT, N_EDGE };
    // An array to allow the edges to be accessed in the correct order.
    static constexpr std::array<Edge, N_EDGE> edges = { BOTTOM, RIGHT, TOP, LEFT };
    std::array<std::string, N_EDGE> edgeNames = { "bottom", "right", "top", "left" };

    MPI_Comm mpiComm;
    int mpiSize = 0;
    int mpiMyRank = -1;
    int localCornerX, localCornerY;
    int localExtentX, localExtentY;
    int globalExtentX, globalExtentY;
    // mpi rank ID and extent for each edge direction
    std::array<std::vector<int>, N_EDGE> neighbourRanks;
    std::array<std::vector<int>, N_EDGE> neighbourExtents;
    std::array<std::vector<int>, N_EDGE> neighbourHaloSend;
    std::array<std::vector<int>, N_EDGE> neighbourHaloRecv;
    std::array<std::vector<int>, N_EDGE> neighbourRanksPeriodic;
    std::array<std::vector<int>, N_EDGE> neighbourExtentsPeriodic;
    std::array<std::vector<int>, N_EDGE> neighbourHaloSendPeriodic;
    std::array<std::vector<int>, N_EDGE> neighbourHaloRecvPeriodic;
#endif

private:
    /*!
     * @brief Read neighbour metadata from netCDF file
     *
     * @param netCDF file with partition metadata
     */
    void readNeighbourData(netCDF::NcFile& ncFile);

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
    bool isCartesian = false;
    // Are the more complex coordinates stored?
    bool hasParameters = false;
#ifdef USE_MPI
    const std::string bboxName = "bounding_boxes";
    const std::string neighbourName = "connectivity";
#endif
};

} /* namespace Nextsim */

#endif /* MODELMETADATA_HPP */
