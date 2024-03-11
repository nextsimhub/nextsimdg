/*!
 * @file IStructure.hpp
 *
 * @date Dec 17, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#ifndef ISTRUCTURE_HPP
#define ISTRUCTURE_HPP

#include "include/ModelMetadata.hpp"
#include "include/ModelState.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <string>

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
    virtual ~IStructure() = default;

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
        const std::string& filePath, ModelMetadata& metadata)
        = 0;
#else
    virtual ModelState getModelState(const std::string& filePath) = 0;
#endif

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
    virtual void dumpModelState(
        const ModelState& state, const ModelMetadata& metadata, const std::string& filePath, bool isRestart) const = 0;

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
};

}
#endif /* ISTRUCTURE_HPP */
