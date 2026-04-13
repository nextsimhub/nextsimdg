/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 */

#include "include/ModelMetadata.hpp"

#include "include/Finalizer.hpp"
#include "include/IStructure.hpp"
#include "include/Logged.hpp"
#include "include/ModelMPI.hpp"
#include "include/NextsimModule.hpp"
#ifdef USE_MPI
#include "include/ParallelNetcdfFile.hpp"
#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif
#include "include/Halo.hpp"
#endif
#include "include/gridNames.hpp"
#include <array>
#include <cstddef>
#include <functional>
#include <vector>

#include <ncDim.h>
#ifdef USE_MPI
#include "mpi.h"
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

        std::array<std::string, 2> cornerSuffixes = { "_neighbour_ids", "_neighbour_send" };
        if (btype == periodic) {
            for (auto& suffix : cornerSuffixes) {
                suffix += "_periodic";
            }
        }

        for (auto corner : corners) {
            size_t nStart = 0; // start index in the netCDF variable
            size_t count = 0; // number of elements to read for this rank
            std::vector<int> numCorners = std::vector<int>(mpiSize, 0);
            std::vector<int> offsets = std::vector<int>(mpiSize, 0);
            std::vector<std::reference_wrapper<std::vector<int>>> arrays;

            // pick the correct set of corner arrays (periodic / non‑periodic)
            if (btype == nonPeriodic) {
                arrays = { cornerRanks[corner], cornerHaloSend[corner] };
            } else { // periodic
                arrays = { cornerRanksPeriodic[corner], cornerHaloSendPeriodic[corner] };
            }

            // variable that stores *how many* corner entries each rank has
            varName = cornerNames[corner] + "_neighbours";
            if (btype == periodic) {
                varName += "_periodic";
            }
            neighbourGroup.getVar(varName).getVar(
                { 0 }, { static_cast<size_t>(mpiSize) }, numCorners.data());

            // compute the global offset for this rank
            MPI_Exscan(&numCorners[mpiMyRank], &nStart, 1, MPI_INT, MPI_SUM, modelMPI.getComm());
            if (mpiMyRank == 0) {
                nStart = 0; // MPI_Exscan undefined on rank 0 → set manually
            }
            count = numCorners[mpiMyRank];

            if (count) {
                // allocate and read each corner‑related array
                for (size_t i = 0; i < arrays.size(); ++i) {
                    arrays[i].get().resize(count, 0);
                    varName = cornerNames[corner] + cornerSuffixes[i];
                    neighbourGroup.getVar(varName).getVar(
                        { nStart }, { count }, arrays[i].get().data());
                }
            }
        }
    }
}

void ModelMetadata::getPartitionMetadata(std::string partitionFile)
{
    try {
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
        if (!globalExtentX) {
            globalExtentX = ncFile.getDim("NX").getSize();
        } else if (globalExtentX != ncFile.getDim("NX").getSize()) {
            throw std::runtime_error("ModelMetadata: Inconsistent global x-extent between "
                                     "partition and input files.");
        }
        if (!globalExtentY) {
            globalExtentY = ncFile.getDim("NY").getSize();
        } else if (globalExtentX != ncFile.getDim("NY").getSize()) {
            throw std::runtime_error("ModelMetadata: Inconsistent global y-extent between "
                                     "partition and input files.");
        }
        netCDF::NcGroup bboxGroup(ncFile.getGroup(bboxName));

        std::vector<size_t> rank(1, modelMPI.getRank());
        bboxGroup.getVar("domain_x").getVar(rank, &localCornerX);
        bboxGroup.getVar("domain_y").getVar(rank, &localCornerY);
        bboxGroup.getVar("domain_extent_x").getVar(rank, &localExtentX);
        bboxGroup.getVar("domain_extent_y").getVar(rank, &localExtentY);

        readNeighbourData(ncFile);

        // cornerHaloRecv doesn't need to be read because it can be easily calculated.
        for (auto corner : corners) {
            if (cornerRanks[corner].size()) {
                cornerHaloRecv[corner].resize(1);
                cornerHaloRecv[corner][0] = 2 * (localExtentX + localExtentY) + corner;
            }
            if (cornerRanksPeriodic[corner].size()) {
                cornerHaloRecvPeriodic[corner].resize(1);
                cornerHaloRecvPeriodic[corner][0] = 2 * (localExtentX + localExtentY) + corner;
            }
        }

        // gather rank extents in X & Y direction for all processes
        rankExtentsX.resize(modelMPI.getSize(), 0);
        rankExtentsY.resize(modelMPI.getSize(), 0);
        MPI_Allgather(
            &localExtentX, 1, MPI_INT, rankExtentsX.data(), 1, MPI_INT, modelMPI.getComm());
        MPI_Allgather(
            &localExtentY, 1, MPI_INT, rankExtentsY.data(), 1, MPI_INT, modelMPI.getComm());

        ncFile.close();
    } catch (netCDF::exceptions::NcException& e) {
        std::cerr << "Failed to open partition file [" << partitionFile << "] :: " << e.what()
                  << std::endl;
        auto& modelMPI = ModelMPI::getInstance();
        MPI_Abort(modelMPI.getComm(), 1);
    }
}

int ModelMetadata::getLocalCornerX() const { return localCornerX; }
int ModelMetadata::getLocalCornerY() const { return localCornerY; }
int ModelMetadata::getLocalExtentX() const { return localExtentX; }
int ModelMetadata::getLocalExtentY() const { return localExtentY; }
int ModelMetadata::getGlobalExtentX() const { return globalExtentX; }
int ModelMetadata::getGlobalExtentY() const { return globalExtentY; }
const std::vector<int>& ModelMetadata::getRankExtentsX() const { return rankExtentsX; }
const std::vector<int>& ModelMetadata::getRankExtentsY() const { return rankExtentsY; }
#else

ModelMetadata::ModelMetadata()
{
    isInitialized = true;
    static bool doneOnce = doOnce();
}

#endif

void ModelMetadata::setDimensionsFromFile(const std::string& filename)
{
    if (filename.empty()) {
        throw std::runtime_error(
            "ModelMetadata :: setDimensionsFromFile() called without input file.");
    }

    // Set to record the names of dimensions found in the file
    std::set<std::string> dimNames;

    try {
#ifdef USE_MPI
        auto& modelMPI = ModelMPI::getInstance();
        netCDF::NcFilePar ncFile(filename, netCDF::NcFile::read, modelMPI.getComm());
#else
        netCDF::NcFile ncFile(filename, netCDF::NcFile::read);
#endif

        // Dimensions and DG components
        std::multimap<std::string, netCDF::NcDim> dimMap = ncFile.getDims();
        for (auto& entry : ModelArray::definedDimensions) {
            const ModelArray::Dimension& dimType = entry.first;
            if (dimType == ModelArray::Dimension::DG || dimType == ModelArray::Dimension::DGSTRESS
                || dimType == ModelArray::Dimension::NCOORDS) {
                // TODO: Assert that DG in the file equals the compile time DG in the model (#205)
                continue;
            }

            const ModelArray::DimensionSpec& dimensionSpec = entry.second;
            // Find dimensions in the netCDF file by their name in the ModelArray details
            netCDF::NcDim dim = ncFile.getDim(dimensionSpec.name);
            // Also check the old name
            if (dim.isNull()) {
                dim = ncFile.getDim(dimensionSpec.altName);
            }
            if (dim.isNull()) {
                Logged::warning(
                    "ModelMetadata: No netCDF dimension found corresponding to the dimension named "
                    + dimensionSpec.name + " or " + dimensionSpec.altName);
                continue;
            } else {
                dimNames.insert(dimensionSpec.name);
            }
#ifdef USE_MPI
            size_t localLength;
            size_t start;
            if (dimType == ModelArray::Dimension::X) {
                localLength = getLocalExtentX();
                start = getLocalCornerX();
            } else if (dimType == ModelArray::Dimension::Y) {
                localLength = getLocalExtentY();
                start = getLocalCornerY();
            } else if (dimType == ModelArray::Dimension::XVERTEX) {
                localLength = getLocalExtentX() + 1;
                start = getLocalCornerX();
            } else if (dimType == ModelArray::Dimension::YVERTEX) {
                localLength = getLocalExtentY() + 1;
                start = getLocalCornerY();
            } else if (dimType == ModelArray::Dimension::XCG) {
                localLength = CGDEGREE * getLocalExtentX() + 1;
                start = CGDEGREE * getLocalCornerX();
            } else if (dimType == ModelArray::Dimension::YCG) {
                localLength = CGDEGREE * getLocalExtentY() + 1;
                start = CGDEGREE * getLocalCornerY();
            } else {
                localLength = dim.getSize();
                start = 0;
            }
            ModelArray::setDimension(
                dimType, dim.getSize(), localLength + 2 * Halo::haloWidth, start);
#else
            ModelArray::setDimension(dimType, dim.getSize());
#endif
        }
    } catch (const netCDF::exceptions::NcException& nce) {
        std::string ncWhat(nce.what());
        ncWhat += ": " + filename;
        throw std::runtime_error(ncWhat);
    }

    // Throw an error if we didn't find any dimensions
    if (dimNames.empty()) {
        throw std::out_of_range(
            "ModelMetadata: No netCDF dimensions in input file '" + filename + "'");
    }
}

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
    xiosHandler.setCalendarStart(time);
#endif
}

void ModelMetadata::incrementTime(const Duration& step)
{
    m_time += step;
#ifdef USE_XIOS
    Xios& xiosHandler = Xios::getInstance();
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
