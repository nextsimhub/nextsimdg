/*!
 * @file IStructure.hpp
 *
 * @date Dec 17, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef CORE_SRC_INCLUDE_ISTRUCTURE_HPP
#define CORE_SRC_INCLUDE_ISTRUCTURE_HPP

#include "include/ElementData.hpp"

#include <boost/algorithm/string/predicate.hpp>
#include <ncFile.h>
#include <ncGroup.h>
#include <string>

// See https://isocpp.org/wiki/faq/pointers-to-members#macro-for-ptr-to-memfn
#define CALL_MEMBER_FN(object, ptrToMember) ((object).*(ptrToMember))

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
    virtual ~IStructure() = default;

    /*!
     * @brief Dumps the data to a file path.
     *
     * @param filePath The path to attempt writing the data to.
     */
    inline void init(const std::string& filePath)
    {
        netCDF::NcFile ncFile(filePath, netCDF::NcFile::FileMode::read);
        init(ncFile);
        ncFile.close();
    }

    /*!
     * @brief Initializes the structure based on the contents of the structure
     * group of the input file.
     *
     * @param grp The NetCDF group instance that holds the structure
     * information.
     */
    virtual void init(const netCDF::NcGroup& grp)
    {
        netCDF::NcGroup metaGroup(grp.getGroup(metadataNodeName));
        netCDF::NcGroup dataGroup(grp.getGroup(dataNodeName));

        initMeta(metaGroup);
        initData(dataGroup);
    }

    /*!
     * @brief Initializes the structure of the IStructure from metadata.
     *
     * @param metaGroup The NetCDF group instance holding the structure
     * metadata.
     */
    virtual void initMeta(const netCDF::NcGroup& metaGroup) = 0;

    /*!
     * @brief Initializes the contents of the IStructure from data.
     *
     * @param dataGroup The NetCDF group instance holding the data.
     */
    virtual void initData(const netCDF::NcGroup& metaGroup) = 0;

    //! Returns the structure name that this class will process
    std::string structureType() const { return processedStructureName; }
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
    virtual void dumpMeta(netCDF::NcGroup& metaGroup) const
    {
        metaGroup.putAtt(typeNodeName, processedStructureName);
    }

    /*!
     * @brief Dumps the prognostic data to a netCDF node.
     *
     * @param dataGroup The top-level node to write the data to.
     */
    virtual void dumpData(netCDF::NcGroup& dataGroup) const = 0;

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
        netCDF::NcGroup metaGroup = headGroup.addGroup(metadataNodeName);
        netCDF::NcGroup dataGroup = headGroup.addGroup(dataNodeName);
        dumpMeta(metaGroup);
        dumpData(dataGroup);
    }

    /*!
     * @brief Dumps the data to a file path.
     *
     * @param filePath The path to attempt writing the data to.
     */
    inline void dump(const std::string& filePath)
    {
        netCDF::NcFile ncFile(filePath, netCDF::NcFile::FileMode::replace);
        dump(ncFile);
        ncFile.close();
    }

    /*!
     * @brief Resets the data cursor.
     *
     * @returns 0 as an int.
     */
    virtual int resetCursor() { return 0; };

    /*!
     * @brief Returns whether the current cursor is valid.
     */
    virtual bool validCursor() const { return false; };

    /*!
     * @brief Returns the data value at the cursor.
     */
    virtual ElementData& cursorData() = 0;

    /*!
     * @brief Returns a const reference to the data value at the cursor.
     */
    virtual const ElementData& cursorData() const = 0;

    /*!
     * @brief Increments the cursor.
     */
    virtual void incrCursor() = 0;

    class Cursor {
    public:
        virtual ~Cursor() = default;
        virtual IStructure& operator=(const int) const = 0;
        virtual operator bool() const = 0;
        virtual ElementData& operator*() const = 0;
        virtual ElementData* operator->() const = 0;
        virtual IStructure& operator++() const = 0;
    };

protected:
    //! Name of the metadata node.
    std::string metadataNodeName = "structure";
    //! Name of the data node.
    std::string dataNodeName = "data";
    //! Name of the structure type processed by this class.
    std::string processedStructureName = "none";
    //! The name of the node holding the name of the structure type processed
    //! by this class.
    const std::string typeNodeName = "type";
};

}
#endif /* CORE_SRC_INCLUDE_ISTRUCTURE_HPP */
