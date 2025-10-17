/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 */

#include "include/ModelMetadata.hpp"

#include "include/Finalizer.hpp"
#include "include/IStructure.hpp"
#include "include/ModelMPI.hpp"
#include "include/NextsimModule.hpp"
#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif
#include "include/gridNames.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

#ifdef USE_MPI
#include "mpi.h"
#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>
#endif
#ifdef USE_OASIS
#include <oasis_c.h>
#endif

namespace Nextsim {

const std::string& ModelMetadata::structureName() const
{
    return Module::getImplementation<IStructure>().structureType();
}

#ifdef USE_MPI
ModelMetadata::ModelMetadata(std::string partitionFile)
{
    if (partitionFile.empty()) {
        throw std::runtime_error(
            "ModelMetadata :: getInstance() called without partition file in MPI build.");
    }
    getPartitionMetadata(partitionFile);
    static bool doneOnce = doOnce();
    isInitialized = true;
}

void ModelMetadata::readNeighbourData(netCDF::NcFile& ncFile)
{
    netCDF::NcGroup neighbourGroup(ncFile.getGroup(neighbourName));
    std::string varName {};
    auto& modelMPI = ModelMPI::getInstance();
    auto mpiSize = modelMPI.getSize();
    auto mpiMyRank = modelMPI.getRank();
    enum BoundaryType { nonPeriodic, periodic };
    for (BoundaryType btype : { nonPeriodic, periodic }) {
        // Use btype as needed
        std::array<std::string, 4> suffixes = { "_neighbour_ids", "_neighbour_halos",
            "_neighbour_halo_send", "_neighbour_halo_recv" };
        if (btype == periodic) {
            for (auto& suffix : suffixes) {
                suffix += "_periodic";
            }
        }
        for (auto edge : edges) {
            size_t nStart = 0; // start point in metadata arrays
            size_t count = 0; // number of elements to read from metadata arrays
            std::vector<int> numNeighbours = std::vector<int>(mpiSize, 0);
            std::vector<int> offsets = std::vector<int>(mpiSize, 0);
            std::vector<std::reference_wrapper<std::vector<int>>> arrays;

            if (btype == nonPeriodic) {
                arrays = { neighbourRanks[edge], neighbourExtents[edge], neighbourHaloSend[edge],
                    neighbourHaloRecv[edge] };
            } else if (btype == periodic) {
                arrays = { neighbourRanksPeriodic[edge], neighbourExtentsPeriodic[edge],
                    neighbourHaloSendPeriodic[edge], neighbourHaloRecvPeriodic[edge] };
            }

            varName = edgeNames[edge] + "_neighbours";
            if (btype == periodic) {
                varName += "_periodic";
            }
            neighbourGroup.getVar(varName).getVar(
                { 0 }, { static_cast<size_t>(mpiSize) }, numNeighbours.data());

            // compute start index for each process
            MPI_Exscan(&numNeighbours[mpiMyRank], &nStart, 1, MPI_INT, MPI_SUM, modelMPI.getComm());
            if (mpiMyRank == 0) {
                // MPI_Exscan is undefined on the first rank. So to be safe we manually set nStart
                // to 0. (see e.g., https://www.open-mpi.org/doc/v4.1/man3/MPI_Exscan.3.php)
                nStart = 0;
            }
            // how many elements to read for each process
            count = numNeighbours[mpiMyRank];

            if (count) {
                // initialize neighbour info to zero
                for (size_t i = 0; i < arrays.size(); ++i) {
                    arrays[i].get().resize(count, 0);
                    varName = edgeNames[edge] + suffixes[i];
                    neighbourGroup.getVar(varName).getVar(
                        { nStart }, { count }, arrays[i].get().data());
                }
            }
        }
    }
}

void ModelMetadata::getPartitionMetadata(std::string partitionFile)
{
    netCDF::NcFile ncFile(partitionFile, netCDF::NcFile::read);
    int sizes = ncFile.getDim("L").getSize();
    int nBoxes = ncFile.getDim("P").getSize();
    auto& modelMPI = ModelMPI::getInstance();
    auto mpiSize = modelMPI.getSize();
    if (nBoxes != mpiSize) {
        std::string errorMsg = "Number of MPI ranks " + std::to_string(mpiSize) + " <> "
            + std::to_string(nBoxes) + "\n";
        throw std::runtime_error(errorMsg);
    }
    globalExtentX = ncFile.getDim("NX").getSize();
    globalExtentY = ncFile.getDim("NY").getSize();
    netCDF::NcGroup bboxGroup(ncFile.getGroup(bboxName));

    std::vector<size_t> rank(1, modelMPI.getRank());
    bboxGroup.getVar("domain_x").getVar(rank, &localCornerX);
    bboxGroup.getVar("domain_y").getVar(rank, &localCornerY);
    bboxGroup.getVar("domain_extent_x").getVar(rank, &localExtentX);
    bboxGroup.getVar("domain_extent_y").getVar(rank, &localExtentY);

    readNeighbourData(ncFile);

    ncFile.close();
}

int ModelMetadata::getLocalCornerX() const { return localCornerX; }
int ModelMetadata::getLocalCornerY() const { return localCornerY; }
int ModelMetadata::getLocalExtentX() const { return localExtentX; }
int ModelMetadata::getLocalExtentY() const { return localExtentY; }
int ModelMetadata::getGlobalExtentX() const { return globalExtentX; }
int ModelMetadata::getGlobalExtentY() const { return globalExtentY; }
#else

ModelMetadata::ModelMetadata()
{
    isInitialized = true;
    static bool doneOnce = doOnce();
}

#endif

const ModelState& ModelMetadata::extractCoordinates(const ModelState& state)
{
    // More sophisticated grids include both vertex coordinates and grid azimuth values.
    if (state.data.count(coordsName) > 0) {
        m_vertexCoords = state.data.at(coordsName);
        m_gridAzimuth = state.data.at(gridAzimuthName);
        hasParameters = true;
    } else {
        // else don't resize the arrays
        hasParameters = false;
    }

    isCartesian = state.data.count(xName);
    if (isCartesian) {
        m_coord1 = state.data.at(xName);
        m_coord2 = state.data.at(yName);
    } else {
        m_coord1 = state.data.at(longitudeName);
        m_coord2 = state.data.at(latitudeName);
    }

    return state;
}

ModelState& ModelMetadata::affixCoordinates(ModelState& state) const
{
    if (hasParameters) {
        state.data[coordsName] = m_vertexCoords;
        state.data[gridAzimuthName] = m_gridAzimuth;
    }

    if (isCartesian) {
        state.data[xName] = m_coord1;
        state.data[yName] = m_coord2;
    } else {
        state.data[longitudeName] = m_coord1;
        state.data[latitudeName] = m_coord2;
    }
    return state;
}

#ifdef USE_OASIS
void ModelMetadata::initOasis(const bool writeOasisGrid)
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
    if (writeOasisGrid) {
        /* This needs to be figured out, but it's not a priority. Grid writing is
         * not necessary for the type of coupling we'll start with.

        const std::string gridName = "nxts";

        int flag = 1;
        OASIS_CHECK_ERR(oasis_c_start_grids_writing(&flag));

        OASIS_CHECK_ERR(oasis_c_write_grid(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, lon, lat, OASISPartitionId));
        OASIS_CHECK_ERR(oasis_c_write_corner(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, clo, cla, OASISPartitionId));
        OASIS_CHECK_ERR(oasis_c_write_area(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, area, OASISPartitionId));
        OASIS_CHECK_ERR(oasis_c_write_mask(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, angle, OASISPartitionId));

        std::string companion = "land area fraction";
        OASIS_CHECK_ERR(oasis_c_write_frac(
                            gridName.c_str(), nx, ny, nx_loc, ny_loc, mask, OASISPartitionId),
            companion.c_str());
        companion = "land sea mask";
        OASIS_CHECK_ERR(oasis_c_write_mask(
                            gridName.c_str(), nx, ny, nx_loc, ny_loc, mask, OASISPartitionId),
            companion.c_str());
        */
    }
}
#endif

void ModelMetadata::setTimes(const TimePoint& start, const TimePoint& stop, const Duration& step)
{
    this->start = start;
    this->stop = stop;
    this->step = step;
    this->run = stop - start;
    setTime(start);
}

void ModelMetadata::setTimes(const TimePoint& start, const Duration& runLen, const Duration& step)
{
    this->start = start;
    this->stop = start + runLen;
    this->step = step;
    this->run = runLen;
    setTime(start);
}

void ModelMetadata::setTime(const TimePoint& time)
{
    m_time = time;
#ifdef USE_XIOS
    Xios& xiosHandler = Xios::getInstance();
    if (!xiosHandler.isInitialized()) {
        throw std::runtime_error("ModelMetadata: Xios handler has not been initialized");
    }
    xiosHandler.setCalendarStart(time);
#endif
}

void ModelMetadata::incrementTime(const Duration& step)
{
    m_time += step;
#ifdef USE_XIOS
    Xios& xiosHandler = Xios::getInstance();
    if (!xiosHandler.isInitialized()) {
        throw std::runtime_error("ModelMetadata: Xios handler has not been initialized");
    }
    xiosHandler.incrementCalendar();
#endif
}

void ModelMetadata::finalize() { }

bool ModelMetadata::doOnce()
{
    // Register the finalization function here
    Finalizer::registerUnique(finalize);
    return true;
}

} /* namespace Nextsim */
