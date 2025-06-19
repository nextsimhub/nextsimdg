/*!
 * @file ModelMetadata.hpp
 *
 * @date 19 Jun 2025
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
private:
#ifdef USE_MPI
    ModelMetadata(std::string partitionFile);
#else
    ModelMetadata();
#endif
    // Prevent copying
    ModelMetadata(const ModelMetadata&) = delete;
    //! Performs some one-time initialization for the class. Returns true.
    static bool doOnce();

public:
#ifdef USE_MPI
    inline static ModelMetadata& getInstance(std::string partitionFile = "")
    {
        static ModelMetadata instance = ModelMetadata(partitionFile);
        if (instance.isInitialized) {
            return instance;
        } else {
            throw std::runtime_error("ModelMetadata :: Object needs to be initialized before use.");
        }
    }
#else
    inline static ModelMetadata& getInstance()
    {
        static ModelMetadata instance = ModelMetadata();
        if (instance.isInitialized) {
            return instance;
        } else {
            throw std::runtime_error("ModelMetadata :: Object needs to be initialized before use.");
        }
    }
#endif

    // finalize ModelMetadata
    static void finalize();

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
    /*!
     * @brief Extracts and sets MPI partition metadata from partition file
     *
     * @param file with partition metadata
     */
    void getPartitionMetadata(std::string partitionFile);

    /*!
     * @brief Gets the local X coordinate of the partition's lower-left corner.
     * @return The X index of the local partition's corner.
     */
    int getLocalCornerX() const;
    /*!
     * @brief Gets the local Y coordinate of the partition's lower-left corner.
     * @return The Y index of the local partition's corner.
     */
    int getLocalCornerY() const;
    /*!
     * @brief Gets the extent of the partition in the X direction.
     * @return The number of grid points in X for the local partition.
     */
    int getLocalExtentX() const;
    /*!
     * @brief Gets the extent of the partition in the Y direction.
     * @return The number of grid points in Y for the local partition.
     */
    int getLocalExtentY() const;
    /*!
     * @brief Gets the global extent of the grid in the X direction.
     * @return The total number of grid points in X for the global domain.
     */
    int getGlobalExtentX() const;
    /*!
     * @brief Gets the global extent of the grid in the Y direction.
     * @return The total number of grid points in Y for the global domain.
     */
    int getGlobalExtentY() const;

    enum Edge { BOTTOM, RIGHT, TOP, LEFT, N_EDGE };
    // An array to allow the edges to be accessed in the correct order.
    static constexpr std::array<Edge, N_EDGE> edges = { BOTTOM, RIGHT, TOP, LEFT };
    std::array<std::string, N_EDGE> edgeNames = { "bottom", "right", "top", "left" };

    typedef std::array<std::vector<int>, N_EDGE> neighbourArray;
    neighbourArray neighbourRanks;
    neighbourArray neighbourExtents;
    neighbourArray neighbourHaloSend;
    neighbourArray neighbourHaloRecv;
    neighbourArray neighbourRanksPeriodic;
    neighbourArray neighbourExtentsPeriodic;
    neighbourArray neighbourHaloSendPeriodic;
    neighbourArray neighbourHaloRecvPeriodic;
#endif

private:
    /*!
     * @brief Read neighbour metadata from netCDF file
     *
     * @param netCDF file with partition metadata
     */
    void readNeighbourData(netCDF::NcFile& ncFile);

#ifdef USE_MPI
    void setMpiMetadata(MPI_Comm comm);
#endif

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
    // has metadata been initialized
    bool isInitialized = false;
#ifdef USE_MPI
    int localCornerX, localCornerY;
    int localExtentX, localExtentY;
    int globalExtentX, globalExtentY;
    const std::string bboxName = "bounding_boxes";
    const std::string neighbourName = "connectivity";
#endif
};

} /* namespace Nextsim */

#endif /* MODELMETADATA_HPP */
