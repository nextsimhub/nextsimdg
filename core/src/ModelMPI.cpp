/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifdef USE_MPI
#include "include/ModelMPI.hpp"

namespace Nextsim {

ModelMPI::ModelMPI(MPI_Comm comm)
    : m_comm(comm)
{
    MPI_Comm_size(m_comm, &m_size);
    MPI_Comm_rank(m_comm, &m_rank);
    static bool doneOnce = doOnce();
    isInitialized = true;
}

MPI_Comm ModelMPI::getComm() const { return m_comm; }

int ModelMPI::getSize() const { return m_size; }

int ModelMPI::getRank() const { return m_rank; }

bool ModelMPI::doOnce()
{
    // Register the finalization function here
    Finalizer::registerUnique(finalize);
    return true;
}

void ModelMPI::finalize() { }

} // namespace Nextsim

#endif // USE_MPI
