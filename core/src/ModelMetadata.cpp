/*!
 * @file ModelMetadata.cpp
 *
 * @date Jun 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelMetadata.hpp"

#include "include/StructureModule.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {

const std::string& ModelMetadata::structureName() const
{
    return Module::getImplementation<IStructure>().structureType();
}

#ifdef USE_MPI
void ModelMetadata::setMpiMetadata(MPI_Comm comm)
{
    mpiComm = comm;
    MPI_Comm_size(mpiComm, &mpiSize);
    MPI_Comm_rank(mpiComm, &mpiMyRank);
}
#endif

const ModelState& ModelMetadata::extractCoordinates(const ModelState& state)
{
    if (state.data.count(coordsName) > 0) {
        m_vertexCoords = state.data.at(coordsName);
    }

    isCartesian = state.data.count(xName);
    if (isCartesian) {
        m_coord1 = state.data.at(xName);
        m_coord2 = state.data.at(yName);
        m_gridAzimuth.resize();
        m_gridAzimuth = 0;
    } else {
        m_coord1 = state.data.at(longitudeName);
        m_coord2 = state.data.at(latitudeName);
        m_gridAzimuth = state.data.at(gridAzimuthName);
    }

    return state;
}

ModelState& ModelMetadata::affixCoordinates(ModelState& state) const
{
    if (m_vertexCoords.trueSize() > 0) {
        state.data[coordsName] = m_vertexCoords;
    }

    if (isCartesian) {
        state.data[xName] = m_coord1;
        state.data[yName] = m_coord2;
    } else {
        state.data[longitudeName] = m_coord1;
        state.data[latitudeName] = m_coord2;
        state.data[gridAzimuthName] = m_gridAzimuth;
    }
    return state;
}
} /* namespace Nextsim */
