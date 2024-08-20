/*!
 * @file ModelMetadata.cpp
 *
 * @date 20 August 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelMetadata.hpp"

#include "include/StructureModule.hpp"
#include "include/gridNames.hpp"

#ifdef USE_MPI
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
ModelMetadata::ModelMetadata(std::string partitionFile, MPI_Comm comm)
{
    setMpiMetadata(comm);
    getPartitionMetadata(partitionFile);
}

void ModelMetadata::setMpiMetadata(MPI_Comm comm)
{
    mpiComm = comm;
    MPI_Comm_size(mpiComm, &mpiSize);
    MPI_Comm_rank(mpiComm, &mpiMyRank);
}

void ModelMetadata::getPartitionMetadata(std::string partitionFile)
{
    // TODO: Move the reading of the partition file to its own class
    netCDF::NcFile ncFile(partitionFile, netCDF::NcFile::read);
    int sizes = ncFile.getDim("L").getSize();
    int nBoxes = ncFile.getDim("P").getSize();
    if (nBoxes != mpiSize) {
        std::string errorMsg = "Number of MPI ranks " + std::to_string(mpiSize) + " <> "
            + std::to_string(nBoxes) + "\n";
        throw std::runtime_error(errorMsg);
    }
    globalExtentX = ncFile.getDim("NX").getSize();
    globalExtentY = ncFile.getDim("NY").getSize();
    netCDF::NcGroup bboxGroup(ncFile.getGroup(bboxName));
    std::vector<size_t> index(1, mpiMyRank);
    bboxGroup.getVar("global_x").getVar(index, &localCornerX);
    bboxGroup.getVar("global_y").getVar(index, &localCornerY);
    bboxGroup.getVar("local_extent_x").getVar(index, &localExtentX);
    bboxGroup.getVar("local_extent_y").getVar(index, &localExtentY);
    ncFile.close();
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
} /* namespace Nextsim */
