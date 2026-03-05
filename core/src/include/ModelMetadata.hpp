/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief Model metadata as a singleton class
 * @details
 *
 * Model metadata as a singleton class
 *
 * This class holds the metadata pertaining to the model as a whole, as well as the MPI metadata
 * which is used for halo exchange.
 *
 * The first time getInstance() is called, the class is initialized. If compiled with MPI then a
 * partition file must be supplied.
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

    /*!
     * @brief   Set dimensions information based on the contents of an input file.
     *
     * @param   filename the name of the file
     * @details If an input file hasn't been read yet, the dimensions are read from the file and
     *          set. Otherwise, a consistency check is made against the dimensions read from file
     *          and already set.
     */
    void setDimensionsFromFile(const std::string& filename);

    // finalize ModelMetadata
    static void finalize();

    /*! Sets the model start, stop and step times directly.
     *
     * @param start. The model initial TimePoint.
     * @param stop. The model final TimePoint.
     * @param step. The model advection/thermodynamics step Duration.
     */
    void setTimes(const TimePoint& start, const TimePoint& stop, const Duration& step);

    /*! Sets the model start, stop and step times from a run length.
     *
     * @param start. The model initial TimePoint.
     * @param runLength. The model run Duration.
     * @param step. The model advection/thermodynamics step Duration.
     */
    void setTimes(const TimePoint& start, const Duration& runLength, const Duration& step);

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

    //! Returns the model start time
    inline const TimePoint& startTime() const { return start; }
    //! Returns the model stop time
    inline const TimePoint& stopTime() const { return stop; }
    //! Returns the model step length
    inline const Duration& stepLength() const { return step; }
    //! Returns the model run length
    inline const Duration& runLength() const { return run; }

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
    /*!
     * @brief Gets the extents of the grid in the X direction for all ranks.
     * @return A vector containing the number of grid points in X for each rank.
     */
    std::vector<int> getRankExtentsX() const;
    /*!
     * @brief Gets the extents of the grid in the Y direction for all ranks.
     * @return A vector containing the number of grid points in Y for each rank.
     */
    std::vector<int> getRankExtentsY() const;

    enum Edge { BOTTOM, RIGHT, TOP, LEFT, N_EDGE };
    // An array to allow the edges to be accessed in the correct order.
    static constexpr std::array<Edge, N_EDGE> edges = { BOTTOM, RIGHT, TOP, LEFT };
    std::array<std::string, N_EDGE> edgeNames = { "bottom", "right", "top", "left" };

    enum Corner { TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT, N_CORNER };
    static constexpr std::array<Corner, N_CORNER> corners
        = { TOP_LEFT, TOP_RIGHT, BOTTOM_RIGHT, BOTTOM_LEFT };
    std::array<std::string, N_CORNER> cornerNames
        = { "top_left", "top_right", "bottom_right", "bottom_left" };

    using NeighbourArray = std::array<std::vector<int>, N_EDGE>;
    NeighbourArray neighbourRanks;
    NeighbourArray neighbourExtents;
    NeighbourArray neighbourHaloSend;
    NeighbourArray neighbourHaloRecv;
    NeighbourArray neighbourRanksPeriodic;
    NeighbourArray neighbourExtentsPeriodic;
    NeighbourArray neighbourHaloSendPeriodic;
    NeighbourArray neighbourHaloRecvPeriodic;

    using CornerArray = std::array<std::vector<int>, N_CORNER>;
    CornerArray cornerRanks;
    CornerArray cornerHaloSend;
    CornerArray cornerHaloRecv;
    CornerArray cornerRanksPeriodic;
    CornerArray cornerHaloSendPeriodic;
    CornerArray cornerHaloRecvPeriodic;
#endif

    std::string initialFileName;
    std::string finalFileName;
    // Period between restart file outputs
    Duration restartPeriod;

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

    TimePoint start;
    TimePoint stop;
    Duration step;
    Duration run;

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
    std::vector<int> rankExtentsX; // vector containing domain extents for each rank x-direction
    std::vector<int> rankExtentsY; // vector containing domain extents for each rank y-direction
    const std::string bboxName = "bounding_boxes";
    const std::string neighbourName = "connectivity";
#endif
};

} /* namespace Nextsim */

#endif /* MODELMETADATA_HPP */
