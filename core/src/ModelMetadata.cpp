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

class CoordError : public std::domain_error {
public:
    CoordError(const std::string& type, const std::string& coord)
    : std::domain_error(std::string("Model is configured as " + type + ", no valid " + coord + " available"))
    {
    }
};
const ModelArray& ModelMetadata::longitude() const
{
    if (isCartesian)
        throw CoordError("Cartesian", "longitude");
    else
        return m_coord1;
}

const ModelArray& ModelMetadata::latitude() const
{
    if (isCartesian)
        throw CoordError("Cartesian", "latitude");
    else
        return m_coord2;
}

const ModelArray& ModelMetadata::x() const
{
    if (!isCartesian)
        throw CoordError("spherical", "x");
    else
        return m_coord1;
}

const ModelArray& ModelMetadata::y() const
{
    if (!isCartesian)
        throw CoordError("spherical", "y");
    else
        return m_coord2;
}
} /* namespace Nextsim */
