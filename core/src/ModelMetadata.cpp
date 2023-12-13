/*!
 * @file ModelMetadata.cpp
 *
 * @date Jun 29, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelMetadata.hpp"

#include "include/StructureModule.hpp"

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

} /* namespace Nextsim */
