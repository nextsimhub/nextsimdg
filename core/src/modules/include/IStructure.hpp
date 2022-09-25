/*!
 * @file IStructure.hpp
 *
 * @date Dec 17, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ISTRUCTURE_HPP
#define ISTRUCTURE_HPP

#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <string>

#ifdef USE_MPI
#include "include/MpiUtils.hpp"
#endif // USE_MPI

namespace Nextsim {

/*!
 * @brief Interface class for the model structure.
 *
 * @details Derived classes will hold the structured data of the model. Member
 * functions also provide
 * · input and output of restart data and data fields for output, either
 *  natively or on a reshaped grid
 * · iteration over the element data, using a given function.
 * This should allow derived classes to implement both Eulerian grids and
 * Lagrangian meshes.
 */
class IStructure {
public:
    IStructure() { }

#ifdef USE_MPI
    // FIXME: This ctor cannot be used with Module infrastructure
    IStructure(MPI_Comm comm)
    {
        mpiComm = comm;
        CHECK_MPI(MPI_Comm_size(comm, &mpiSize));
        CHECK_MPI(MPI_Comm_rank(comm, &mpiRank));
    }
#endif // USE_MPI

    virtual ~IStructure() = default;

#ifdef USE_MPI
    /*!
     * @brief Sets the MPI Communicator for this distributed structure.
     */
    void setComm(MPI_Comm comm)
    {
        mpiComm = comm;
        CHECK_MPI(MPI_Comm_size(comm, &mpiSize));
        CHECK_MPI(MPI_Comm_rank(comm, &mpiRank));
    }

    /*!
     * @brief Returns MPI Communicator associated with this structure.
     */
    MPI_Comm getComm() const { return mpiComm; }
#endif // USE_MPI

    /*!
     * @brief Dumps the data to a file path.
     *
     * @param filePath The path to attempt writing the data to.
     */
    //    virtual void init(const std::string& filePath) = 0;

    /*!
     * @brief Returns the ModelState stored in the file
     */
#ifdef USE_MPI
    virtual ModelState getModelState(
        const std::string& restartFilePath, const std::string& partitionFilePath)
        = 0;
#else
    virtual ModelState getModelState(const std::string& filePath) = 0;
#endif // USE_MPI

    //! Returns the structure name that this class will process
    virtual const std::string& structureType() const { return processedStructureName; }
    /*!
     * @brief Checks if the passed string matches (ignoring case) the name of
     * the structure that this class constructs.
     *
     * @param str The string to be checked.
     */
    inline bool structureTypeCheck(const std::string& str) const
    {
        return boost::algorithm::iequals(structureType(), str);
    }

    //! The number of ice layers in this data structure.
    virtual int nIceLayers() const = 0;

    /*!
     * @brief Dumps the data to a file path.
     *
     * @param filePath The path to attempt writing the data to.
     */
    //    virtual void dump(const std::string& filePath) const = 0;

    /*!
     * @brief Dumps the given ModelState to the given file path.
     *
     * @param state The ModelState data
     * @param filePath The path to attempt to write the data to.
     */
    virtual void dumpModelState(const ModelState& state, const ModelMetadata& metadata,
        const std::string& filePath, bool isRestart) const
        = 0;

    // Node names in the default structure

    //! Returns the name of the metadata node.
    static const std::string metadataNodeName() { return "metadata"; }
    //! Returns the name of the data node.
    static const std::string dataNodeName() { return "data"; }
    //! The name of the group holding the definitive structure type
    static const std::string structureNodeName() { return "structure"; }
    //! The name of the node holding the name of the structure type processed
    //! by this class.
    static const std::string typeNodeName() { return "type"; }

private:
    //! Name of the structure type processed by this class.
    const std::string processedStructureName = "none";
#ifdef USE_MPI
    MPI_Comm mpiComm; // MPI communicator
    int mpiSize = 0; // Number of MPI processes in communicator
    int mpiRank = -1; // MPI rank
#endif // USE_MPI
};

}
#endif /* ISTRUCTURE_HPP */
