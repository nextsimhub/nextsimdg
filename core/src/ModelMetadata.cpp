/*!
 * @file ModelMetadata.cpp
 *
 * @date 21 August 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelMetadata.hpp"

#include "include/IStructure.hpp"
#include "include/NextsimModule.hpp"
#include "include/gridNames.hpp"

#ifdef USE_MPI
#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>
#endif

namespace Nextsim {

ModelMetadata::ModelMetadata()
    : m_vertexCoords(ModelArray::Type::VERTEX)
    , m_coord1(ModelArray::Type::H)
    , m_coord2(ModelArray::Type::H)
    , m_gridAzimuth(ModelArray::Type::H)
    , isCartesian(false)
    , hasParameters(false)
    {
    }

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
    bboxGroup.getVar("domain_x").getVar(index, &localCornerX);
    bboxGroup.getVar("domain_y").getVar(index, &localCornerY);
    bboxGroup.getVar("domain_extent_x").getVar(index, &localExtentX);
    bboxGroup.getVar("domain_extent_y").getVar(index, &localExtentY);
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
        state.data.emplace(coordsName, m_vertexCoords);
        state.data.emplace(gridAzimuthName, m_gridAzimuth);
    }

    std::pair<std::string, std::string> coordNames;
    if (isCartesian) {
        coordNames = {xName, yName};
    } else {
        coordNames = {longitudeName, latitudeName};
    }
    state.data.emplace(coordNames.first, m_coord1);
    state.data.emplace(coordNames.second, m_coord2);

    return state;
}
} /* namespace Nextsim */
