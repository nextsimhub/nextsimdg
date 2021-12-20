/*!
 * @file IStructure.hpp
 *
 * @date Dec 17, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_ISTRUCTURE_HPP_
#define CORE_SRC_INCLUDE_ISTRUCTURE_HPP_

#include <boost/algorithm/string/predicate.hpp>
#include <ncGroup.h>
#include "/opt/home/include/ncGroup.h"
#include <string>
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
    virtual ~IStructure();

    /*!
     * @brief Initializes the structure based on the contents of the structure
     * group of the input file.
     *
     * @param grp The NetCDF group instance that holds the structure
     * information.
     */
    virtual void init(netCDF::NcGroup& grp) = 0;

    //! Returns the structure name that this class will process
    virtual std::string structureType() const { return "none"; }
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

    /*!
     * @brief Dumps the structural metadata to a netCDF node.
     *
     * @param metaGroup The top-level node to write the metadata to.
     */
    virtual void dumpMeta(netCDF::NcGroup& metaGroup) const = 0;

    /*!
     * @brief Dumps the prognostic data to a netCDF node.
     *
     * @param dataGroup The top-level node to write the data to.
     */
    virtual void dumpData(netCDF::NcGroup& dataGroup) const = 0;

    /*!
     * @brief Dumps the data and metadata to two netCDF nodes.
     *
     * @param dataGroup The top-level node to write the data to.
     * @param metaGroup The top-level node to write the metadata to.
     */
    inline void dump(netCDF::NcGroup& metaGroup, netCDF::NcGroup& dataGroup) const
    {
        dumpMeta(metaGroup);
        dumpData(dataGroup);
    }

    /*!
     * @brief Dumps the data and metadata to two sub-groups.
     *
     * @details The structure metadata will be dumped to immediate subnodes
     * with names specified by the ::metadataNodeName and ::dataNodeName public
     * member variables.
     *
     * @param headGroup The top-level node to hold the metadata and data nodes.
     */
    inline void dump(netCDF::NcGroup& headGroup) const
    {
        dump(headGroup.addGroup(metadataNodeName), headGroup.addGroup(dataNodeName));
    }

    /*!
     * @brief Dumps the data to a file path.
     *
     * @param filePath The path to attempt writing the data to.
     */
    inline void dump(const std::string& filePath) {
            dump(netCDF::NcFile(filePath, netCDF::NcFile::FileMode::replace));
    }
    //! Default name of the metadata node
    std::string metadataNodeName = "structure";
    std::string dataNodeName = "data";
};

#endif /* CORE_SRC_INCLUDE_ISTRUCTURE_HPP_ */
