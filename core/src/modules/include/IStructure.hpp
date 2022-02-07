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
    IStructure()
        : cursor(*this)
    {
    }
    virtual ~IStructure() = default;

    /*!
     * @brief Dumps the data to a file path.
     *
     * @param filePath The path to attempt writing the data to.
     */
    virtual void init(const std::string& filePath) = 0;

    //! Returns the structure name that this class will process
    virtual std::string structureType() const { return processedStructureName; }
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
    virtual void dump(const std::string& filePath) const = 0;
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
        Cursor(IStructure& ownerer)
            : owner(ownerer)
        {
        }
        ~Cursor() = default;
        IStructure& operator=(const int i) const
        {
            if (0 == i)
                owner.resetCursor();
            return owner;
        }
        operator bool() const { return owner.validCursor(); };
        ElementData& operator*() const { return owner.cursorData(); };
        ElementData* operator->() const { return &owner.cursorData(); };
        IStructure& operator++() const
        {
            owner.incrCursor();
            return owner;
        };

    private:
        IStructure& owner;
    };

    const Cursor cursor;

    // Node names in the default structure

    //! Returns the name of the metadata node.
    static const std::string metadataNodeName() { return "structure"; };
    //! Returns the name of the data node.
    static const std::string dataNodeName() { return "data"; };
    //! The name of the node holding the name of the structure type processed
    //! by this class.
    static const std::string typeNodeName() { return "type"; };

private:
    //! Name of the structure type processed by this class.
    const std::string processedStructureName = "none";
};

}
#endif /* CORE_SRC_INCLUDE_ISTRUCTURE_HPP */
