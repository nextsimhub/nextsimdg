/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief Model MPI singleton
 * @details
 *
 * Model MPI singleton
 *
 * MPI metadata is stored in the singleton.
 */

#ifndef MODELMPI_HPP
#define MODELMPI_HPP
#ifdef USE_MPI

#include "include/Finalizer.hpp"
#include <mpi.h>
#include <stdexcept>

namespace Nextsim {

class ModelMPI {
private:
    ModelMPI(MPI_Comm comm);
    // Prevent copying
    ModelMPI(const ModelMPI&) = delete;
    //! Performs some one-time initialization for the class. Returns true.
    static bool doOnce();

public:
    inline static ModelMPI& getInstance(MPI_Comm comm = MPI_COMM_NULL)
    {
        static ModelMPI instance = ModelMPI(comm);
        if (instance.isInitialized) {
            return instance;
        } else {
            throw std::runtime_error("ModelMPI :: Object needs to be initialized before use.");
        }
    }

    // finalize ModelMPI
    static void finalize();

    /*!
     * @brief Gets the MPI communicator.
     * @return The MPI communicator for all processes.
     */
    MPI_Comm getComm() const;

    /*!
     * @brief Gets the MPI size (number of processes).
     * @return The size of the communicator.
     */
    int getSize() const;

    /*!
     * @brief Gets the MPI rank (process ID).
     * @return The rank of this process in the communicator.
     */
    int getRank() const;

private:
    bool isInitialized = false;
    MPI_Comm m_comm;
    int m_size = 0;
    int m_rank = 0;
};

} /* namespace Nextsim */

#endif /* MPI_HPP */
#endif /* MODELMPI_HPP */
